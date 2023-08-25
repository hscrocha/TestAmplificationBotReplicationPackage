import sys
import csv
import json
import re
import shutil
import zipfile
import glob
import os
import requests
from dateutil import parser as dup
import pandas as pd
import datetime
import hashlib
from collections import Counter

from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
from matplotlib_venn import venn3, venn3_circles, venn3_unweighted
from matplotlib import pyplot as plt

sys.path.insert(1, '../../small-amp/runner/')
import reports as rp

"""
How to:

First save a valid github token to ~/.smallAmpCI
Then, in this folder run a python interactive
python3

Then import the py script
import smallampCI

Then load the artifacts:
smallampCI.loadArtifacts('mabdi', 'Roassal3', 1712957423, 1)
Args -> git user, project, run id, attempt
It will download the artifacts and save in the ../experiments and temp/ folders.

Use this command to extract the zips:
python3 smallampCI.py tmp/smallAmp-overview-XXX.zip

Use this command to generate csv data:
python3 smallampCI.py tmp/XXX > tmp/XXX.csv

Use this command to copy the csv text:
cat tmp/XXX.csv | pbcopy

Then go to a google spreadsheet and paste the data

"""

runs = [
      'Seaside-run8',
      'Seaside-run10',
      'Seaside-run11',
      'Seaside-run9',
      'Seaside-run12',
      'Seaside-run13',
      'PolyMath-run5',
      'PolyMath-run7',
      'PolyMath-run8',
      'PolyMath-run6',
      'PolyMath-run9',
      'PolyMath-run10',
      'NovaStelo-run6',
      'NovaStelo-run8',
      'NovaStelo-run9',
      'NovaStelo-run7',
      'NovaStelo-run11',
      'NovaStelo-run12',
      'Moose-run4',
      'Moose-run6',
      'Moose-run7',
      'Moose-run5',
      'Moose-run8',
      'Moose-run9',
      'zinc-run5',
      'zinc-run7',
      'zinc-run8',
      'zinc-run6',
      'zinc-run9',
      'zinc-run10'
   ]

def load_workflow_params(run):
   with open(os.path.join('tmp/{}'.format(run), 'workflow_params.txt'), "r") as f:
        workflow_params = json.loads(f.read())
   return      workflow_params

def load_df_from_csv(csvName, alist):
   dtype=   {
            # 'Name': str,
            # 'Grade': int
             'className':str, 'numberOfTestMethodsBeforeShreding':pd.Int64Dtype(),
       'targetLoc':pd.Int64Dtype(), 'mutationScoreBefore':float, 'n_notCoveredInOriginal':pd.Int64Dtype(),
       'numberOfAllMutationsInOriginal':pd.Int64Dtype(), 'testClassTimeToRunInMillis':pd.Int64Dtype(), 'mode':str,
       'folder':str, 'status':str, 'numberOfOriginalTestMethods':pd.Int64Dtype(),
        'n_amplifiedMethods':pd.Int64Dtype(), 'n_newCovered':pd.Int64Dtype(),
       'mutationScoreAfter':float, 'mutationImprove':float, 'timeBudgetFinished': bool,
       'numberOfProcessedMethods':pd.Int64Dtype(), 'duplicateMutants':pd.Int64Dtype(), 'timeTotal':pd.Int64Dtype()
            }
   na_values = ['None', 'NA']

   frames = [ ]

   for run in alist:
      csv_name = 'csvs/{}/{}.csv'.format(run, csvName)
      
      workflow_params = load_workflow_params(run) 
   
      df = pd.read_csv(csv_name, dtype=dtype, na_values=na_values, on_bad_lines= 'warn')
      projectName = run.split('-')[0]
      df['projectName'] =  projectName
      df['mode'] =  workflow_params['mode']
      df['run'] =  run
      frames.append( df )

   df = pd.concat(frames)
   # df_amp = df.where(df['status'] == 'Finished successfully')
   # df_amp = df_amp[df_amp['projectName'].notnull() ]
   df['timeTotal_delta'] = pd.to_timedelta(df['timeTotal'], unit='sec')
   df['originalTestClass'] = df['className'].replace("Test\d+To\d+", "Test", regex=True)
   # df['timeBudgetFinished'] = df['timeBudgetFinished'].apply(lambda x: not pd.isna(x) and x.strip().lower() == 'true')
   return df


def load_error_df_from_csv(csvName, alist):
   dtype=   {
         'className': str, 'Error': str, 'Info': str
      }
   na_values = ['None', 'NA']

   frames = [ ]

   for run in alist:
      csv_name = 'csvs/{}/{}.csv'.format(run, csvName)
      with open(os.path.join('tmp/{}'.format(run), 'workflow_params.txt'), "r") as f:
        workflow_params = json.loads(f.read())
   
      df = pd.read_csv(csv_name, dtype=dtype, na_values=na_values, on_bad_lines= 'warn')
      projectName = run.split('-')[0]
      df['projectName'] =  projectName
      df['originalTestClass'] = df['className'].replace("Test\d+To\d+", "Test", regex=True)
      df['mode'] =  workflow_params['mode']
      df['run'] =  run
      frames.append( df )

   df = pd.concat(frames)
   return df

def succ_fixed_df():
   return load_df_from_csv('succ_fixed', runs)

def succ_all_df():
   return load_df_from_csv('succ_not_fixed', runs)


def error_all_df():
   return load_error_df_from_csv('fail_not_fixed', runs)


def r_all():
   df = succ_all_df()
   return df[ df['mode'] == 'fseRank' ]

def r_fixed():
   df = succ_fixed_df()
   return df[ df['mode'] == 'fseRank' ]


def rq1_forRun(errors, success, success_fixed):
   
   return {
      'n_classes_no_shred': len( set( errors['originalTestClass'] ).union( set( success_fixed['originalTestClass'] ) ) ),
      'n_classes_after_shred': len(errors ) + len( success  ),
      'n_testmethods': success_fixed['numberOfOriginalTestMethods'].sum(),
      'n_mutants': success_fixed['numberOfAllMutationsInOriginal'].sum(),
      'n_alive_original': success_fixed['n_notCoveredInOriginal'].sum(),
      'n_smallamp_run': len (success) + len(errors),
      'n_smallamp_finished': len (success),
      'n_SANoUncovered': len(errors[errors['Error'] == 'SANoUncovered']),
      'n_SANoGreen': len(errors[errors['Error'] == 'SANoGreenTest']),
      'n_SAFlaky': len(errors[errors['Error'] == 'SAFlakyMutationTesting']),
      'n_unknown_err': len(errors[errors['Error'] == 'Unknown']),
      'n_improved': len(success[success['n_newCovered'] > 0] ),
      'percent_improved': 100 * len(success[success['n_newCovered'] > 0] )/ len(success),
      'n_amplifiedMethods': success['n_amplifiedMethods'].sum(),
      'n_newCovered': success_fixed['n_newCovered'].sum(),
      'increased_kill': 100* success_fixed['n_newCovered'].sum() / (success_fixed['numberOfAllMutationsInOriginal'] - success_fixed['n_notCoveredInOriginal']).sum(),
      'workflow_time': 'TODO',
      'n_timeouts': len (success[success['timeBudgetFinished'] ] ),
      'n_skipped_methods': (success[success['timeBudgetFinished'] ]['numberOfOriginalTestMethods'] - 
            success[success['timeBudgetFinished'] ]['numberOfProcessedMethods'] ).sum(),
      'n_methodsin_timeout': success[success['timeBudgetFinished'] ]['numberOfOriginalTestMethods'],
      'n_large_classes': len( success_fixed[success_fixed['numberOfTestMethodsBeforeShreding'] > 15 ])
   }


def workflow_run_duration_sum(run_list):
   x = datetime.timedelta(0)
   for item in run_list:
      x = x + workflow_run_duration(item)      
   return  x

def workflow_run_duration(run):
   raw_data_folder = '../experiments'
   project = run.split('-')[0]
   
   for r_id in os.listdir(os.path.join(raw_data_folder, project)):
      if os.path.exists(os.path.join(raw_data_folder, project,r_id, 'run_status.json')):
         with open(os.path.join(raw_data_folder, project,r_id, 'run_status.json'), "r") as f:
            run_status_json = json.loads(f.read())
         if str(run_status_json['run_number']) == run.split('-')[1][len('run'):]:
            return dup.parse(run_status_json['updated_at']) - dup.parse(run_status_json['run_started_at'])
   print('failed to find duration for ' + run)
   return '?'

def count_crashes_list(run_list):
   v = 0
   for run in run_list:
      v = v +  count_crashes(run)
   return v

def count_crashes(run):
   raw_data_folder = '../experiments'
   text_to_search = b'Deleting heartbeat file.'
   project = run.split('-')[0]
   n_recovers = 0
   for r_id in os.listdir(os.path.join(raw_data_folder, project)):
      if os.path.exists(os.path.join(raw_data_folder, project,r_id, 'run_status.json')):
         with open(os.path.join(raw_data_folder, project,r_id, 'run_status.json'), "r") as f:
            run_status_json = json.loads(f.read())
         if str(run_status_json['run_number']) == run.split('-')[1][len('run'):]:
            z = os.path.join(raw_data_folder, project,r_id, 'run_logs.zip')
            if os.path.exists(z):
               unzipped_file = zipfile.ZipFile(z, "r")
               for job_id in range(0,8):   
                  a_file = unzipped_file.read("Job number {}/4_Run mabdismallamp-action@main.txt".format(job_id))
                  # print(run,,a_file.count(text_to_search))
                  n_recovers += a_file.count(text_to_search)
   return n_recovers         


def flatten3(ll):
   return [t for ttt in ll for tt in ttt for t in tt]

def find_textin_logs(run,text):
   raw_data_folder = '../experiments'
   text_to_search = text
   project = run.split('-')[0]
   found_cases = []
   for r_id in os.listdir(os.path.join(raw_data_folder, project)):
      if os.path.exists(os.path.join(raw_data_folder, project,r_id, 'run_status.json')):
         with open(os.path.join(raw_data_folder, project,r_id, 'run_status.json'), "r") as f:
            run_status_json = json.loads(f.read())
         if str(run_status_json['run_number']) == run.split('-')[1][len('run'):]:
            z = os.path.join(raw_data_folder, project,r_id, 'run_logs.zip')
            if os.path.exists(z):
               unzipped_file = zipfile.ZipFile(z, "r")
               for job_id in range(0,8):   
                  a_file = unzipped_file.read("Job number {}/4_Run mabdismallamp-action@main.txt".format(job_id))
                  # print(run,,a_file.count(text_to_search))
                  # n_recovers += a_file.count(text_to_search)
                  found = re.findall(text, a_file.decode("UTF-8"), flags=re.MULTILINE)
                  if found:
                     found_cases.append( found )
   return found_cases         



def rq2():
   run_groups = [ runs[i::6] for i in range(0, 6)]
   errors = error_all_df()
   succ_fixed = succ_fixed_df()
   succ_all = succ_all_df()

   print("""\\begin{table*}
\\caption{The result of the quantitative experiment.}
\\label{table:results}

\\begin{tabular}{ |p{5cm}|p{1.5cm}p{1.5cm}p{1.5cm}|p{1.5cm}p{1.5cm}p{1.5cm}|} 
 \\hline
 & \\multicolumn{3}{p{4.5cm}|}{\\textbf{Prioritizing Enabled}} & \\multicolumn{3}{p{4.5cm}|}{\\textbf{No Prioritizing}} \\\\
 & \\textbf{\\#1} & \\textbf{\\#2}& \\textbf{\\#3} & \\textbf{\\#1} & \\textbf{\\#2} & \\textbf{\\#3} \\\\
 \\hline""")


  
   tmp = [ len( set( errors[errors['run'].isin(g)]['originalTestClass'] ).union( set( set( succ_all[succ_all['run'].isin(g)]['originalTestClass'] ) ) )) for g in run_groups]
   print('\\# Class before sharding & {} & {} & {} & {} & {} & {}  \\\\'.format( * tmp ))
   tmp = [ len( set( errors[errors['run'].isin(g)]['className'] ).union( set( set( succ_all[succ_all['run'].isin(g)]['className'] ) ) )) for g in run_groups]
   print('\\# Class after sharding & {} & {} & {} & {} & {} & {}  \\\\'.format( * tmp  ))
   tmp = [ len(errors[errors['run'].isin(g) & (errors['Error'] == 'SANoUncovered')])  for g in run_groups]
   print('\\# No Uncovered & {} & {} & {} & {} & {} & {}  \\\\'.format(* tmp  ))
   tmp = [ len(errors[errors['run'].isin(g) & (errors['Error'] == 'SANoGreenTest')])  for g in run_groups]
   print('\\# No Green & {} & {} & {} & {} & {} & {}  \\\\'.format(* tmp  ))
   tmp = [ len(errors[errors['run'].isin(g) & (errors['Error'] == 'SAFlakyMutationTesting')])  for g in run_groups]
   print('\\# Flaky MutationTesting& {} & {} & {} & {} & {} & {}  \\\\'.format(* tmp  ))
   tmp = [ len(errors[errors['run'].isin(g) & (errors['Error'] == 'Unknown')])  for g in run_groups]
   print('\\# Unknown error & {} & {} & {} & {} & {} & {}  \\\\'.format(* tmp  ))
   tmp = [ len(succ_all[succ_all['run'].isin(g)])  for g in run_groups]
   print('\\# Executions finished & {} & {} & {} & {} & {} & {}  \\\\'.format(* tmp  ))
   tmp = [ succ_fixed[succ_fixed['run'].isin(g)]['numberOfOriginalTestMethods'].sum()  for g in run_groups]
   print('\\# Test methods & {} & {} & {} & {} & {} & {}  \\\\'.format( * tmp  ))
   tmp = [ succ_fixed[succ_fixed['run'].isin(g)]['numberOfAllMutationsInOriginal'].sum()  for g in run_groups]
   print('\\# Mutants & {} & {} & {} & {} & {} & {}  \\\\'.format( * tmp  ))
   tmp = [ succ_fixed[succ_fixed['run'].isin(g)]['n_notCoveredInOriginal'].sum()  for g in run_groups]
   print('\\# Mutants alive original & {} & {} & {} & {} & {} & {}  \\\\'.format(* tmp  ))
   tmp = [ len(succ_all[succ_all['run'].isin(g) & (succ_all['n_newCovered'] > 0)])  for g in run_groups]
   print('\\# Executions having improvement & {} & {} & {} & {} & {} & {}  \\\\'.format(* tmp  ))
   tmp = [ 100 * len(succ_all[succ_all['run'].isin(g) & (succ_all['n_newCovered'] > 0)]) / len(succ_all[succ_all['run'].isin(g)])  for g in run_groups]
   print('\\% Executions having improvement & {:.2f}\\% & {:.2f}\\% & {:.2f}\\% & {:.2f}\\% & {:.2f}\\% & {:.2f}\\%  \\\\'.format(* tmp  ))
   tmp = [ succ_all[succ_all['run'].isin(g) & (succ_all['n_newCovered'] > 0)]['n_amplifiedMethods'].sum() for g in run_groups]
   print('\\# Generated tests & {} & {} & {} & {} & {} & {}  \\\\'.format(* tmp  ))
   tmp = [ succ_all[succ_all['run'].isin(g)]['n_newCovered'].sum() for g in run_groups]
   print('\\# Newly killed mutants & {} & {} & {} & {} & {} & {}  \\\\'.format(* tmp  ))
   tmp = [ 100 * succ_all[succ_all['run'].isin(g)]['n_newCovered'].sum() / (succ_all[succ_all['run'].isin(g)]['numberOfAllMutationsInOriginal'] - succ_all[succ_all['run'].isin(g)]['n_notCoveredInOriginal']).sum() for g in run_groups]
   print('\\% Increased kills & {:.2f}\\% & {:.2f}\\% & {:.2f}\\% & {:.2f}\\% & {:.2f}\\% & {:.2f}\\% \\\\'.format(* tmp  ))
   tmp = [ count_crashes_list(g) for g in run_groups]
   print('\\# Recovered crashes & {} & {} & {} & {} & {} & {}  \\\\'.format(* tmp  ))
   tmp = [ len(succ_fixed[succ_fixed['run'].isin(g) & ( succ_fixed['numberOfOriginalTestMethods'] > 15)])  for g in run_groups]
   print('\\# Large test classes & {} & {} & {} & {} & {} & {}  \\\\'.format(* tmp  ))
   tmp = [ len(succ_all[succ_all['run'].isin(g) & ( succ_all['timeBudgetFinished'])]) for g in run_groups]
   print('\\# Time budget finished & {} & {} & {} & {} & {} & {}  \\\\'.format(* tmp  ))
   tmp = [ (succ_all[succ_all['run'].isin(g)]['numberOfOriginalTestMethods'] - succ_all[succ_all['run'].isin(g)]['numberOfProcessedMethods']).sum() for g in run_groups]
   print('\\# Test methods skipped & {} & {} & {} & {} & {} & {}  \\\\'.format(* tmp  ))
   
   tmp = [ workflow_run_duration_sum(g) for g in run_groups]
   print('Workflow duration: All & {} & {} & {} & {} & {} & {}  \\\\'.format(* tmp  ))
   

   for runItem in run_groups[0]:
      projectName = runItem.split('-')[0]
      
      tmp = [ workflow_run_duration(r) for g in run_groups for r in g if r.startswith(projectName)]
      if projectName == 'zinc':
         projectName = 'Zinc'
      tmp.insert(0, projectName)
      print('\\hspace{{1em}}{} & {} & {} & {} & {} & {} & {}  \\\\'.format(* tmp  ))
   
      
   
   print("""  \\hline
\\end{tabular}
\\end{table*}""")


def rq3_new():
   run_groups = [ runs[i::6] for i in range(0, 6)]
   errors = error_all_df()
   succ_fixed = succ_fixed_df()
   succ_all = succ_all_df()

   rgg = run_groups[0]

   projects = [ runItem.split('-')[0] for runItem in run_groups[0] ]
   print("""\\begin{table}
\\caption{The state of test classes in the dataset.}
\\label{table:tclasses}
\\begin{small}
\\begin{tabular}{ p{4cm}|p{0.5cm}p{0.5cm}p{0.5cm}p{0.5cm}p{0.5cm}} 
 %\\hline
 """)
   
   print("& \\rotatebox{{30}}{{\\textbf{{{}}}}} & \\rotatebox{{30}}{{\\textbf{{{}}}}}& \\rotatebox{{30}}{{\\textbf{{{}}}}}& \\rotatebox{{30}}{{\\textbf{{{}}}}}& \\rotatebox{{30}}{{\\textbf{{{}}}}} \\\\\\hline".format(* projects))



   r1 = [ len( set( errors[errors['run'] == r]['originalTestClass'] ).union( set( set( succ_all[succ_all['run']==r]['originalTestClass'] ) ) )) for projectName in projects for r in rgg if r.startswith(projectName)]
   print('\\# Classes before sharding& {} & {} & {} & {} & {} \\\\'.format(* r1  ))
   r2 = [ len( set( errors[errors['run'].eq(r) & (errors['className'] != errors['originalTestClass'])]['originalTestClass'] ).union( set( set( succ_fixed[succ_fixed['run'].eq(r) & ( succ_fixed['className'] != succ_fixed['originalTestClass'])]['originalTestClass'])))) for projectName in projects for r in rgg if r.startswith(projectName)]
   print('\\# Large test classes & {} & {} & {} & {} & {} \\\\'.format(* r2  ))
   r3 = [ len( set( errors[errors['run'] == r]['className'] ).union( set( set( succ_all[succ_all['run'] == r]['className'] ) ) )) for projectName in projects for r in rgg if r.startswith(projectName)]
   print('\\# Classes after sharding& {} & {} & {} & {} & {} \\\\'.format(* r3  ))
   r4 = [ len(errors[(errors['run']==r) & (errors['Error'] == 'SANoUncovered')]) for projectName in projects for r in rgg if r.startswith(projectName)]
   print('\\# Classes with 100\% coverage & {} & {} & {} & {} & {}  \\\\'.format(* r4  ))
   r5 = [ len(errors[(errors['run']==r) & (errors['Error'] == 'SANoGreenTest')]) for projectName in projects for r in rgg if r.startswith(projectName)]
   print('\\# Classes without any green test & {} & {} & {} & {} & {}  \\\\'.format(* r5  ))

   # r6 = [ len( set(errors[errors['run'].eq(r) & (errors['Error'] == 'SAFlakyMutationTesting')]['originalTestClass']) )  for projectName in projects for r in rgg if r.startswith(projectName)]

   print('\\hline')
   r_min = [sum(r) for r in zip(r4,r5) ] 
   r_sum = [r[0] - r[1] for r in zip(r3,r_min) ]
   print('\\# Classes to be amplified & {} & {} & {} & {} & {}  \\\\'.format(* r_sum  ))
   
   
   print("""  %\\hline
\\end{tabular}
\\end{small}
\\end{table}""")


   print("""





      """)



   print("""\\begin{table*}
\\caption{The result of the quantitative experiment.}
\\label{table:results}

\\begin{tabular}{ |p{5cm}|p{1.5cm}p{1.5cm}p{1.5cm}|p{1.5cm}p{1.5cm}p{1.5cm}|} 
 \\hline
 & \\multicolumn{3}{p{4.5cm}|}{\\textbf{Prioritizing Enabled}} & \\multicolumn{3}{p{4.5cm}|}{\\textbf{No Prioritizing}} \\\\
 & \\textbf{\\#1} & \\textbf{\\#2}& \\textbf{\\#3} & \\textbf{\\#1} & \\textbf{\\#2} & \\textbf{\\#3} \\\\
 \\hline""")
   tmp = [ len(errors[errors['run'].isin(g) & (errors['Error'] == 'SAFlakyMutationTesting')])  for g in run_groups]
   print('\\# Flaky MutationTesting& {} & {} & {} & {} & {} & {}  \\\\'.format(* tmp  ))
   tmp = [ len(errors[errors['run'].isin(g) & (errors['Error'] == 'Unknown')])  for g in run_groups]
   print('\\# Unknown error & {} & {} & {} & {} & {} & {}  \\\\'.format(* tmp  ))
   tmp = [ len(succ_all[succ_all['run'].isin(g)])  for g in run_groups]
   print('\\# Executions finished & {} & {} & {} & {} & {} & {}  \\\\'.format(* tmp  ))
   tmp = [ succ_fixed[succ_fixed['run'].isin(g)]['numberOfOriginalTestMethods'].sum()  for g in run_groups]
   print('\\# Test methods & {} & {} & {} & {} & {} & {}  \\\\'.format( * tmp  ))
   tmp = [ succ_fixed[succ_fixed['run'].isin(g)]['numberOfAllMutationsInOriginal'].sum()  for g in run_groups]
   print('%\\# Mutants & {} & {} & {} & {} & {} & {}  \\\\'.format( * tmp  ))
   tmp = [ succ_fixed[succ_fixed['run'].isin(g)]['n_notCoveredInOriginal'].sum()  for g in run_groups]
   print('%\\# Mutants alive original & {} & {} & {} & {} & {} & {}  \\\\'.format(* tmp  ))
   tmp = [ len(succ_all[succ_all['run'].isin(g) & (succ_all['n_newCovered'] > 0)])  for g in run_groups]
   print('\\# Executions having improvement & {} & {} & {} & {} & {} & {}  \\\\'.format(* tmp  ))
   tmp = [ 100 * len(succ_all[succ_all['run'].isin(g) & (succ_all['n_newCovered'] > 0)]) / len(succ_all[succ_all['run'].isin(g)])  for g in run_groups]
   print('\\% Executions having improvement & {:.2f}\\% & {:.2f}\\% & {:.2f}\\% & {:.2f}\\% & {:.2f}\\% & {:.2f}\\%  \\\\'.format(* tmp  ))
   tmp = [ succ_all[succ_all['run'].isin(g) & (succ_all['n_newCovered'] > 0)]['n_amplifiedMethods'].sum() for g in run_groups]
   print('\\# Generated tests & {} & {} & {} & {} & {} & {}  \\\\'.format(* tmp  ))
   tmp = [ succ_all[succ_all['run'].isin(g)]['n_newCovered'].sum() for g in run_groups]
   print('%\\# Newly killed mutants & {} & {} & {} & {} & {} & {}  \\\\'.format(* tmp  ))
   tmp = [ 100 * succ_all[succ_all['run'].isin(g)]['n_newCovered'].sum() / (succ_all[succ_all['run'].isin(g)]['numberOfAllMutationsInOriginal'] - succ_all[succ_all['run'].isin(g)]['n_notCoveredInOriginal']).sum() for g in run_groups]
   print('%\\% Increased kills & {:.2f}\\% & {:.2f}\\% & {:.2f}\\% & {:.2f}\\% & {:.2f}\\% & {:.2f}\\% \\\\'.format(* tmp  ))
   tmp = [ count_crashes_list(g) for g in run_groups]
   print('\\# Recovered crashes & {} & {} & {} & {} & {} & {}  \\\\'.format(* tmp  ))
   tmp = [ len(succ_all[succ_all['run'].isin(g) & ( succ_all['timeBudgetFinished'])]) for g in run_groups]
   print('\\# Time budget finished & {} & {} & {} & {} & {} & {}  \\\\'.format(* tmp  ))
   tmp = [ (succ_all[succ_all['run'].isin(g)]['numberOfOriginalTestMethods'] - succ_all[succ_all['run'].isin(g)]['numberOfProcessedMethods']).sum() for g in run_groups]
   print('\\# Test methods skipped & {} & {} & {} & {} & {} & {}  \\\\'.format(* tmp  ))
   
   tmp = [ workflow_run_duration_sum(g) for g in run_groups]
   print('Workflow duration: All & {} & {} & {} & {} & {} & {}  \\\\'.format(* tmp  ))
   

   for runItem in run_groups[0]:
      projectName = runItem.split('-')[0]
      
      tmp = [ workflow_run_duration(r) for g in run_groups for r in g if r.startswith(projectName)]
      if projectName == 'zinc':
         projectName = 'Zinc'
      tmp.insert(0, projectName)
      print('\\hspace{{1em}}{} & {} & {} & {} & {} & {} & {}  \\\\'.format(* tmp  ))
   
      
   
   print("""  \\hline
\\end{tabular}
\\end{table*}""")


def rq_crashes():
   run_groups = [ runs[i::6] for i in range(0, 6)]
   errors = error_all_df()
   succ_fixed = succ_fixed_df()
   succ_all = succ_all_df()


def rqX():
   #TODO create x3 columns per run + avg
   succ_fixed = succ_fixed_df()
   succ_fixed = succ_fixed[ succ_fixed['mode'] == 'fseRank' ]
   errors = error_all_df()
   errors = errors[ errors['mode'] == 'fseRank' ]
   succ_all = succ_all_df()
   succ_all = succ_all[ succ_all['mode'] == 'fseRank' ]
   
   first_run_path = runs[0::6]
   first_errors = errors[errors['run'].isin(first_run_path)]
   first_success_fixed = succ_fixed[ succ_fixed['run'].isin(first_run_path) ]
   first_success = succ_all[ succ_all['run'].isin(first_run_path) ]
   


   print('# class before shred\t\t{} X3'.format( len( set( first_errors['originalTestClass'] ).union( set( first_success_fixed['originalTestClass'] ) ) ) ))
   print('# class after shred\t\t{} X3'.format( len(first_errors ) + len( first_success  ) ) )
   print('# test methods\t\t{} X3'.format( first_success_fixed['numberOfOriginalTestMethods'].sum() ))
   print('# mutants\t\t{} X3'.format( first_success_fixed['numberOfAllMutationsInOriginal'].sum() ))

   print('# mutants alive orig\t\t{} X3'.format( first_success_fixed['n_notCoveredInOriginal'].sum() ))

   print('# executions \t\t{}'.format( len (succ_all) + len(errors) ))
   print('# executions finished\t\t{}'.format(len(succ_all)))
   print('# executions SANoUncovered\t\t{}'.format( len(errors[errors['Error'] == 'SANoUncovered']) ))
   print('# executions SANoGreen\t\t{}'.format( len(errors[errors['Error'] == 'SANoGreenTest'])))
   print('# executions SAFlaky\t\t{}'.format( len(errors[errors['Error'] == 'SAFlakyMutationTesting'])))
   print('# executions Unknown\t\t{}'.format( len(errors[errors['Error'] == 'Unknown'])))

   print('# finished & imporved\t\t{}'.format(len(succ_all[succ_all['n_newCovered'] > 0] )))
   print('% finished & imporved\t\t{:.2f}'.format( 100 * len(succ_all[succ_all['n_newCovered'] > 0] ) / len(succ_all)  ))

   print('# generated tests\t\t{}'.format( first_success['n_amplifiedMethods'].sum() ))
   print('# newly killed mutants\t\t{}'.format( succ_fixed['n_newCovered'].sum() ))
   print('% increased kills\t\t{:.2f}'.format( 100 * succ_fixed['n_newCovered'].sum() / 
            (succ_fixed['numberOfAllMutationsInOriginal'] - succ_fixed['n_notCoveredInOriginal']).sum()) )
   print('workflow time\t\t{}'.format( 'TODO.' ))

   print('# timeouts \t\t{} out of {}'.format( len (succ_all[succ_all['timeBudgetFinished'] ] ),  len(succ_all)))
   print('# test methods skipped because timeouts \t\t{} out of {}'.format(  (succ_all['numberOfOriginalTestMethods'] - succ_all['numberOfProcessedMethods'] ).sum(),  (succ_all['numberOfOriginalTestMethods']).sum()))

   print('# XXL classes\t\t{} out of {}'.format( len( first_success_fixed[first_success_fixed['numberOfTestMethodsBeforeShreding'] > 15 ] ), len(first_success_fixed)))


def rq4():
   run123 = [ runs[0::6], runs[1::6], runs[2::6] ]
   
   for runItem in run123:
      print('=======')
      succ_fixed = succ_fixed_df()
      succ_fixed = succ_fixed[ (succ_fixed['mode'] == 'fseRank') & (succ_fixed['run'].isin(runItem)) ]
      succ_all = succ_all_df()
      succ_all = succ_all[ (succ_all['mode'] == 'fseRank') & (succ_all['run'].isin(runItem))  ]


      shreds = succ_all[succ_all['numberOfTestMethodsBeforeShreding'] > 15]
      shredded_classes = shreds['originalTestClass'].unique()
      duplicates = succ_fixed[ succ_fixed['className'].isin(shredded_classes) ]['duplicateMutants'].sum()
      allMutantsInShreds = shreds['n_newCovered'].sum()
      print('# duplicated mutants in shreds\t\t{} out of {} ({:.2f}%)'.format(  duplicates,allMutantsInShreds, 100*duplicates/allMutantsInShreds  ) )  
   print('=======All=======')
   succ_fixed = succ_fixed_df()
   succ_fixed = succ_fixed[ (succ_fixed['mode'] == 'fseRank') ]
   succ_all = succ_all_df()
   succ_all = succ_all[ (succ_all['mode'] == 'fseRank')  ]


   shreds = succ_all[succ_all['numberOfTestMethodsBeforeShreding'] > 15]
   shredded_classes = shreds['originalTestClass'].unique()
   duplicates = succ_fixed[ succ_fixed['className'].isin(shredded_classes) ]['duplicateMutants'].sum()
   allMutantsInShreds = shreds['n_newCovered'].sum()
   print('# duplicated mutants in shreds\t\t{} out of {} ({:.2f}%)'.format(  duplicates,allMutantsInShreds, 100*duplicates/allMutantsInShreds  ) )  

def load_jsons_in_run(run,filter):
   jsons = []
   for r in run:
      path = 'tmp/{}'.format(r)
      project = r.split('-')[0]
      files = os.listdir(path)

      with open(os.path.join(path,'todo.txt'), "r") as f:
         classNames = f.read().split('\n')

      for className in classNames:
         if len(filter) > 0 and className not in filter:
            continue
         if os.path.exists(os.path.join(path,className+ ".json")):
            with open(os.path.join(path,className+ ".json"), "r") as f:
               jsons.append(json.loads(f.read()))
   return jsons


def mut_to_string(mut):
   strMut = '-'.join(
         [mut['operatorDescription'], mut['class'], mut['operatorClass'], mut['method'], str(mut['mutationStart']), str(mut['mutationEnd']), mut['mutatedNode'] ]
      )
   result = hashlib.sha256(strMut.encode())
   return result.hexdigest()


def newCoveredInRun(run,filter):
   js = load_jsons_in_run(run,filter)
   newCovered = set([])
   for j in js:
      mutants = [mut_to_string(m) for m in j['newCovered']]
      newCovered.update(mutants)
   return newCovered

def report_set_status(newCovered1, newCovered2, newCovered3):
   covered12 = newCovered1.intersection(newCovered2)
   covered13 = newCovered1.intersection(newCovered3)
   covered23 = newCovered2.intersection(newCovered3)

   allCovered = newCovered1.union(newCovered2).union(newCovered3)
   covered123 = newCovered1.intersection(newCovered2).intersection(newCovered3)

   only1 = newCovered1 - newCovered2 - newCovered3
   only2 = newCovered2 - newCovered1 - newCovered3
   only3 = newCovered3 - newCovered2 - newCovered1

   inMoreThanOne = (covered12 -  covered123).union( covered13 -  covered123 ).union( covered23 -  covered123 )


   print('Covered in run1: \t\t{}'.format( len(newCovered1) ))
   print('Covered in run2: \t\t{}'.format( len(newCovered2) ))
   print('Covered in run3: \t\t{}'.format( len(newCovered3) ))

   print('Covered in run1 only: \t\t{}'.format( len(only1) ))
   print('Covered in run2 only: \t\t{}'.format( len(only2) ))
   print('Covered in run3 only: \t\t{}'.format( len(only3) ))

   print('Covered in all (intersect): \t\t{}'.format(  len(covered123) ))
   print('Covered in more than 1: \t\t{}'.format(  len(inMoreThanOne) ))

   print('All covered (u): \t\t{}'.format(  len(allCovered) ))

   v = venn3_unweighted(subsets = (len(newCovered1 - newCovered2 - newCovered3), len(newCovered2  - newCovered1 - newCovered3), len(covered12 -  covered123), 
         len(newCovered3 - newCovered2 - newCovered1), len(covered13 - covered123), len(covered23 - covered123), len(covered123)), 
         set_labels = ('Run 1', 'Run 2', 'Run 3'),  alpha = 0.5);
   # v.get_patch_by_id('111').set_color('orange')
   # plt.figure(figsize=(4,4))
   plt.show()




def rq3():
   run1 = runs[0::6] 
   run2 = runs[1::6]
   run3 = runs[2::6]

   newCovered1 = newCoveredInRun(run1, [])
   newCovered2 = newCoveredInRun(run2, [])
   newCovered3 = newCoveredInRun(run3, [])

   print('Statues in all classes: ')
   report_set_status(newCovered1, newCovered2, newCovered3)

   succ_all = succ_all_df()
   succ_all = succ_all[ succ_all['mode'] == 'fseRank' ]

   timeout1 = succ_all[succ_all['timeBudgetFinished'] & (succ_all['run'].isin(run1) ) ]
   timeout2 = succ_all[succ_all['timeBudgetFinished'] & (succ_all['run'].isin(run2) ) ]
   timeout3 = succ_all[succ_all['timeBudgetFinished'] & (succ_all['run'].isin(run3) ) ]

   timeout1 = set(timeout1['originalTestClass'])
   timeout2 = set(timeout2['originalTestClass'])
   timeout3 = set(timeout3['originalTestClass'])

   timeout_all = timeout1.intersection(timeout2).intersection(timeout3)
   newCoveredTimeOut1 = newCoveredInRun(run1, timeout_all)
   newCoveredTimeOut2 = newCoveredInRun(run2, timeout_all)
   newCoveredTimeOut3 = newCoveredInRun(run3, timeout_all)

   print('# classes timeouted in all \t\t{}'.format( len (timeout_all )))
   print('Statues in timeouts: ')
   report_set_status(newCoveredTimeOut1, newCoveredTimeOut2, newCoveredTimeOut3)




def rq1():
   succ_all = succ_all_df()
   rnk = succ_all[ succ_all['mode'] == 'fseRank' ]
   non = succ_all[ succ_all['mode'] == 'fseNone' ] 

   rnk_timeout = set(rnk[rnk['timeBudgetFinished']]['originalTestClass'])
   non_timeout = set(non[non['timeBudgetFinished']]['originalTestClass'])
   timeouts = []
   for c in rnk_timeout.intersection(non_timeout):
      if len(set(rnk[(rnk['originalTestClass'] == c) & (rnk['timeBudgetFinished']) ]['run'])) == 3 and len(set(non[(non['originalTestClass'] == c) & (non['timeBudgetFinished']) ]['run'])) == 3:
            timeouts.append(c)

   rnk_fin = set(rnk[~rnk['timeBudgetFinished']]['originalTestClass'])
   non_fin = set(non[~non['timeBudgetFinished']]['originalTestClass'])
   
   # rnk_timeout = set(rnk[rnk['timeBudgetFinished']]['className'])
   # non_timeout = set(non[non['timeBudgetFinished']]['className'])

   
   # timeouts = [x for x in timeouts if x.endswith('Test')]
   timeouts1 = rnk_timeout.union(non_timeout)

   fin = rnk_fin.intersection(non_fin) - timeouts1

   trr = rnk[ rnk['originalTestClass'].isin(timeouts) ]
   tnn = non[ non['originalTestClass'].isin(timeouts) ]

   rr = rnk[ rnk['originalTestClass'].isin(fin) ]
   nn = non[ non['originalTestClass'].isin(fin) ]


   
   print('classes timeouted in both runs \t\t{}'.format( [list(rnk[rnk['originalTestClass'] == t]['projectName'])[0] + ' >> ' + t for t in timeouts] ))
   print('# classes timeouted in both runs \t\t{}'.format( len(timeouts )))
   print('# classes timeouted in at least a run \t\t{}'.format( len(timeouts1 )))
   print('# classes not timedout in both runs \t\t{}'.format( len(fin) ) )

   # n_newCovered    sum
   # 'mutationImprove'] >0

   print('Total:')
   print('# new mutants not ranked, not timeout\t\t{}'.format( nn['n_newCovered'].sum() ))
   print('# new mutants ranked, not timeout\t\t{}'.format( rr['n_newCovered'].sum() ))
   print('# new mutants not ranked, timeout\t\t{}'.format( tnn['n_newCovered'].sum() ))
   print('# new mutants ranked, timeout\t\t{}'.format( trr['n_newCovered'].sum() ))

   # print('# success amp not ranked, not timeout\t\t{}'.format( len(nn[nn['mutationImprove'] >0]) ))
   # print('# success amp ranked, not timeout\t\t{}'.format( len(rr[rr['mutationImprove'] >0]) ))
   # print('# success amp not ranked, timeout\t\t{}'.format( len(tnn[tnn['mutationImprove'] >0]) ))
   # print('# success amp ranked, timeout\t\t{}'.format( len(trr[trr['mutationImprove'] >0]) ))

   # print('Per project')



def rq_crashes_recovered():
   run_groups = [ runs[i::6] for i in range(0, 6)]

   errors = error_all_df()
   tmp = errors[errors['Error'].isin(['Unknown', 'SAFlakyMutationTesting']) ].groupby(['Error', 'projectName', 'originalTestClass']).size()
   
   folders = list(errors[errors['Error'].isin(['Unknown']) ]['run'].unique())
   folders = ['tmp/{0}'.format(f) for f in folders]
   json_files = [f + '/' + pos_json for f in folders for pos_json in os.listdir(f) if pos_json.endswith('.json') and pos_json.startswith('crash_evidence_')]
   mutation_testing = []
   elsewhere = []
   for f in json_files:
      with open(os.path.join(f), "r") as ff:
        jsonf = json.loads(ff.read())
      if 'mutant' in jsonf:
         mutation_testing.append(f)
      else:
         elsewhere.append(f)
   print(tmp)      
   print('All crashes: {}, mutation_testing: {}, elsewhere: {}'.format(len(json_files) , len(mutation_testing), len(elsewhere)))


   timeouts = {str(rg):0 for rg in run_groups}
   # timeouts = {g:find_textin_logs(g, 'Amplification Terminated because timeout') for rg in run_groups for g in rg}
   for rg in run_groups:
      for g in rg:
         timeouts[str(rg)] += len(find_textin_logs(g, 'Amplification Terminated because timeout'))
   # print(timeouts)
   tmp = [ timeouts[str(g)] for g in run_groups]
   print('number of timeouts: {} {} {} {} {} {}'.format(* tmp) )

   timeoutFinal = [ find_textin_logs(g, 'final steps\r?\n.{28} Amplification Terminated because timeout') for rg in run_groups for g in rg]
   print('number of timeouts in final: {}'.format(len(flatten3(timeoutFinal))))
   
   # notTimeoutFinal = [ find_textin_logs(g, '(?!final steps)\r?\n.{28} Amplification Terminated because timeout') for rg in run_groups for g in rg]
   # print('number of not timeouts in final: {}'.format(len(flatten3(notTimeoutFinal))))
   # print(notTimeoutFinal)

   crash_counts = [ find_textin_logs(g, 'Number of crashed for this class up to now: \d+') for rg in run_groups for g in rg]
   crash_counts = flatten3(crash_counts)
   print('number of crash recovery: {}'.format(len(crash_counts)))

   print (dict( (l, crash_counts.count(l) ) for l in crash_counts))


   print('-------Crashes from logs:')

   with open('../crashes.txt') as f:
         logcrashs_all = f.readlines()

   logcrashs = [x for x in logcrashs_all if x != '=====New File======']
   cdet = [i for i in range(0, len(logcrashs)) if 'Number of crashed for this class up to now:' in logcrashs[i]]
   print('All crashes recovered: {}'.format(len(cdet)))
   cdet = [logcrashs[i] for i in range(0, len(logcrashs)) if 'Number of crashed for this class up to now: 1' in logcrashs[i]]
   print('All crashes recovered #1: {}'.format(len(cdet)))
   cdet = [logcrashs[i] for i in range(0, len(logcrashs)) if 'Number of crashed for this class up to now: 2' in logcrashs[i]]
   print('All crashes recovered #2: {}'.format(len(cdet)))

   
   cdet = [logcrashs[i] for i in range(0, len(logcrashs)) if i>1 and 'Number of crashed for this class up to now:' in logcrashs[i] and 'because timeout' in logcrashs[i-1] ]
   print('All timeouts recovered: {}'.format(len(cdet)))

   cdet = [logcrashs[i-1][78:-1] for i in range(0, len(logcrashs)) if i>1 and 'Number of crashed for this class up to now:' in logcrashs[i] and 'because timeout' not in logcrashs[i-1] ]
   print('All non-timeouts recovered: {}'.format(len(cdet)))   

   cdet = [logcrashs[i-2][78:-1] for i in range(0, len(logcrashs)) if i>1 and 'Number of crashed for this class up to now:' in logcrashs[i] and 'because timeout' not in logcrashs[i-1] and 'Segmentation' in logcrashs[i-2] ]
   print('All Segmentation fault: {}'.format(len(cdet)))      

   cdet = [logcrashs[i-1].split(':')[4].strip() for i in range(0, len(logcrashs)) if i>1 and 'Number of crashed for this class up to now:' in logcrashs[i] and 'possible crash for className:' in logcrashs[i-1] ]
   print('All classes having a crash: ' + str(Counter(cdet)))
   
   jsons_dict =  {}
   for f in json_files:
      with open(os.path.join(f), "r") as ff:
        jsonf = json.loads(ff.read())
      for clss in set(cdet):
         if clss not in jsons_dict:
            jsons_dict[clss] = []
         if 'testClass' in jsonf and clss in jsonf['testClass']:
            jsons_dict[clss].append(f)
         
   print('All json files: ' + str(jsons_dict))
   # WAContinuationTest -> WAContinuationTest, timeout
   # 
         

def overall_df():
   df = normalized_df()
   print('{} -> {}'.format('total', str(len(df))))
   for group, frame in df.groupby('status'):
      print('{} -> {}'.format(group, str(len(frame))))
   df2 = df[ df['status'] == 'Finished successfully' ] 
   df2_rank= df2[ (df2['mode'] == 'fseRank') ]
   df2_none= df2[ (df2['mode'] == 'fseNone') ]
   n_improved_classes_rank= len(df2_rank[ df2_rank['mutationImprove'] >0 ])
   n_improved_classes_none= len(df2_none[ df2_none['mutationImprove'] >0 ])

   n_improved_classes= len(df2[ df2['mutationImprove'] >0 ])
   print('# has improved: {}'.format( n_improved_classes ))
   print('% has improved: {:.2f}'.format( n_improved_classes/ len(df2) ))
   print('# has improved ranked: {}'.format( n_improved_classes_rank ))
   print('# has improved none: {}'.format( n_improved_classes_none ))
   print('% has improved ranked: {:.2f}'.format( n_improved_classes_rank / len(df2_rank) ))
   print('% has improved none: {:.2f}'.format( n_improved_classes_none / len(df2_none) ))
            

def per_project_df():
   df = normalized_df()
   print('{} -> {}'.format('total', str(len(df))))
   for group1, frame1 in df.groupby('projectName'):
      print(group1)
      print('  Status:')
      for group2, frame2 in frame1.groupby('status'):
         print('    {} -> {}'.format(group2, str(len(frame2))))
      df2 = frame1[ frame1['status'] == 'Finished successfully' ]
      print('  Improvements:')
      n_improved_classes= len(df2[ df2['mutationImprove'] >0 ])
      n_improved_classes_rank= len(df2[ (df2['mode'] == 'fseRank') & (df2['mutationImprove'] >0) ])
      n_improved_classes_none= len(df2[ (df2['mode'] == 'fseNone') & (df2['mutationImprove'] >0) ])
      print('     # has improved: {}'.format( n_improved_classes ))
      print('     % has improved: {:.2f}'.format( n_improved_classes/ len(df2) ))
      print('     # has improved ranked: {}'.format( n_improved_classes_rank ))
      print('     # has improved none: {}'.format( n_improved_classes_none ))
      

def rq1_repeat():
   loadFrom = '../timeout'
   repeated_runs = [
      'Seaside-run14',
      'Seaside-run15',
      'Seaside-run16',
      'Seaside-run17',
      'Seaside-run18',
      'Seaside-run19',
      'PolyMath-run12',
      'PolyMath-run13',
      'PolyMath-run14',
      'PolyMath-run15',
      'PolyMath-run16',
      'PolyMath-run17',
      'NovaStelo-run13',
      'NovaStelo-run14',
      'NovaStelo-run15',
      'NovaStelo-run16',
      'NovaStelo-run17',
      'NovaStelo-run18',
      'zinc-run11',
      'zinc-run12',
      'zinc-run13',
      'zinc-run14',
      'zinc-run15',
      'zinc-run16'
   ]
   df = load_df_from_csv('succ_fixed', repeated_runs)

   
home = os.path.expanduser("~")
tmpDir = 'tmp'
ssaveTo = '../experiments'
token_file = home + '/.smallAmpCI'

# goto folder
# python3
# import smallampCI
# smallampCI.loadArtifacts('mabdi', 'Roassal3', 1712957423, 1)
def loadArtifacts(user, project, runId, attempt, saveTo = ssaveTo):
   try:
      with open(token_file,'r') as f:
         api_token = f.read().strip()
         print('token loaded successfully, len= ' + str(len(api_token)))
   except OSError as err: 
      print("OS error: {0}".format(err))
      sys.exit('create the file '+ token_file+' including a valid Github token')
   
   session = requests.Session()
   session.headers['Authorization'] = 'token %s' % api_token
   
   url_list_artifacts = 'https://api.github.com/repos/{owner}/{repo}/actions/runs/{run_id}/artifacts'.format(owner=user, repo=project, run_id=runId) # get runId from url
   r = session.get(url_list_artifacts)
   artifacts_json = r.json()
   os.makedirs('{saveTo}/{project}/{runId}'.format(saveTo=saveTo, project=project, runId = runId), exist_ok=True)
   sampleZipName = ''
   if 'artifacts' not in artifacts_json:
      print(artifacts_json)
   for artifact in artifacts_json['artifacts']:
      r= session.get(artifact['archive_download_url'], allow_redirects=True)
      filename = '{saveTo}/{project}/{runId}/{filename}.zip'.format(saveTo=saveTo, project=project, runId = runId, filename=artifact['name'])
      open(filename, 'wb').write(r.content)
      print('Artifact saved: '+ filename)
      sampleZipName = filename
   
   runLogs = 'https://api.github.com/repos/{owner}/{repo}/actions/runs/{run_id}/attempts/{attempt_number}/logs'.format(owner=user, repo=project, run_id=runId, attempt_number= attempt)
   r = session.get(runLogs)
   filename = '{saveTo}/{project}/{runId}/run_logs.zip'.format(saveTo=saveTo, project=project, runId = runId)
   open(filename, 'wb').write(r.content)
   print('Logs saved: '+ filename)


   runStatus = 'https://api.github.com/repos/{owner}/{repo}/actions/runs/{run_id}'.format(owner=user, repo=project, run_id=runId)
   r = session.get(runStatus)
   filename = '{saveTo}/{project}/{runId}/run_status.json'.format(saveTo=saveTo, project=project, runId = runId)
   open(filename, 'wb').write(r.content)
   run_status_json = r.json()
   r = session.get(run_status_json['workflow_url'])
   workflow_basename = os.path.basename(r.json()['path'])

   workflow_file = 'https://raw.githubusercontent.com/{owner}/{repo}/{head_sha}/{path}'.format(owner=user, repo=project, head_sha=run_status_json['head_sha'], path=r.json()['path'])
   r = session.get(workflow_file)

   filename = '{saveTo}/{project}/{runId}/{base}'.format(saveTo=saveTo, project=project, runId = runId, base=workflow_basename)
   open(filename, 'wb').write(r.content)
   print('Workflow saved: '+ filename)
   # the length of run should be: r.json()['updated_at'] - r.json()['run_started_at']
   exec_length = dup.parse(run_status_json['updated_at']) - dup.parse(run_status_json['run_started_at'])
   print('Run is finished in %s' % str(exec_length))
   print('Copying zip files to temp')
   os.system('cp {saveTo}/{project}/{runId}/smallAmp-*.zip {destination}/'.format(saveTo=saveTo, project=project, runId = runId, destination = tmpDir))
   print('Run the script with this path: '+ tmpDir + '/' + sampleZipName + '.zip')


def write_csv_fail(data, directory, filename):
   with open(os.path.join(directory, filename), 'w', newline='') as csvfile:
      fieldnames = ['className', 'Error' , 'Info']
      writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

      writer.writeheader()
      for row in data:
         if row['stat'] == 'error':
            writer.writerow({
               'className': row['className'],
               'Error': row['errDet'], 
               'Info': row['lastMethod']
               })
         elif row['stat'] == 'fail':
            writer.writerow({
               'className': row['className'],
               'Error': 'Unknown', 
               'Info': ''
               })

def write_csv_suc(data, directory, filename):
   with open(os.path.join(directory, filename), 'w', newline='') as csvfile:
      fieldnames = [
               'className', 
               'amplifiedClass',
               'targets',
               'mutationScoreBefore',
               'mutationScoreAfter',
               'mutationImprove',
               'numberOfOriginalTestMethods',
               'targetLoc',
               'testLoc',
               'testAmpLoc',
               'n_amplifiedMethods',
               'n_notCoveredInOriginal',
               'n_newCovered',
               'n_notCoveredInAmplified',
               'n_methodsNotProfiled',
               'timeTotal',
               'max_number_of_changes',
               'numberOfProcessedMethods',
               'testClassTimeToRunInMillis',
               'numberOfAllMutationsInOriginal',
               'numberOfTestMethodsBeforeShreding',
               'timeBudgetFinished',
               'duplicateMutants'
      ]
      writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

      writer.writeheader()
      for row in data:
             jsonObj = row['jsonObj']
             targets = ' '.join(jsonObj['targetClasses'])
             
             writer.writerow({
               'className': row['className'],
               'amplifiedClass' : jsonObj['amplifiedClass'],
               'targets': targets,
               'mutationScoreBefore':      "{:.2f}".format(jsonObj['mutationScoreBefore']),
               'mutationScoreAfter':      "{:.2f}".format(jsonObj['mutationScoreAfter']),
               'mutationImprove':      "{:.2f}".format(jsonObj['mutationScoreAfter'] - jsonObj['mutationScoreBefore']),
               'numberOfOriginalTestMethods':      jsonObj.get('numberOfOriginalTestMethods','NA'),
               'targetLoc': jsonObj['targetLoc'],
               'testLoc':     jsonObj['testLoc'],
               'testAmpLoc':      jsonObj['testAmpLoc'],
               'n_amplifiedMethods': len(jsonObj['amplifiedMethods']),
               'n_notCoveredInOriginal': len(jsonObj['notCoveredInOriginal']),
               'n_newCovered': len(jsonObj['newCovered']),
               'n_notCoveredInAmplified': len(jsonObj['notCoveredInAmplified']),
               'n_methodsNotProfiled':      len(jsonObj['methodsNotProfiled']),
               'timeTotal':      jsonObj['timeTotal'],
               'max_number_of_changes':      max(rp.number_of_changes(jsonObj['amplifiedMethods']) or [0]),
               'numberOfProcessedMethods':      jsonObj.get('numberOfProcessedMethods',0),
               'testClassTimeToRunInMillis':      jsonObj.get('testClassTimeToRunInMillis',0),
               'numberOfAllMutationsInOriginal':      jsonObj.get('numberOfAllMutationsInOriginal',0),
               'numberOfTestMethodsBeforeShreding':      jsonObj.get('numberOfTestMethodsBeforeShreding',0),
               'timeBudgetFinished':      jsonObj.get('timeBudgetFinished',False),
               'duplicateMutants':      jsonObj.get('duplicateMutants',0)
            })
         
def create_csvs(directory):
   fixed = rp.reportAmp_backend(directory, True)
   not_fixed = rp.reportAmp_backend(directory, False)


   toFolder = 'csvs/{}'.format(os.path.basename(directory))
   if not os.path.exists(toFolder):
      os.mkdir(toFolder)

   write_csv_suc([row for row in fixed if row['stat'] == 'success' ], toFolder, 'succ_fixed.csv')
   write_csv_suc([row for row in not_fixed if row['stat'] == 'success' ], toFolder, 'succ_not_fixed.csv')
   
   write_csv_fail([row for row in fixed if row['stat'] != 'success' ], toFolder, 'fail_fixed.csv')
   write_csv_fail([row for row in not_fixed if row['stat'] != 'success' ], toFolder, 'fail_not_fixed.csv')
      

def create_all_csvs(alist):
   for item in alist:
      folder = 'tmp/{}'.format(item)
      create_csvs( folder ) 

def process_folder(folder):
   #print("Processing the folder="+ folder)
   headers_row = 'projectName, className, status, amplifiedClass, targets, mutationScoreBefore, mutationScoreAfter, mutationImprove , numberOfOriginalTestMethods, targetLoc, '\
      'testLoc, testAmpLoc, assertionDensityOriginal, assertionDensityAmplified, originalCoverageStatementes,amplifiedCoverageStatementes, originalCoverageBranches,'\
      'amplifiedCoverageBranches, originalCoverageMethods, amplifiedCoverageMethods, n_amplifiedMethods, n_notCoveredInOriginal, n_newCovered, n_notCoveredInAmplified, n_methodsNotProfiled,'\
      'timeTotal, targetChurn, testChurn, directTestingOriginal, number_of_changes, numberOfProcessedMethods, testClassTimeToRunInMillis, n_mutantsInOriginal, numberOfTestMethodsBeforeShreding,'\
      'timeBudgetFinished, duplicateMutants'
   fieldnames=[ x.strip() for x in headers_row.split(',')]
   #print(fieldnames)
   
   with open(os.path.join(folder, 'workflow_params.txt'), "r") as f:
        workflow_params = json.loads(f.read())
   with open(os.path.join(folder, 'overview-amp.txt'), "r") as f:
        reader = csv.DictReader(f, fieldnames=fieldnames)
        amp_overview = list(reader)
   # with open(os.path.join(folder, 'overview-amp-fixed.txt'), "r") as f:
   #      reader = csv.DictReader(f, fieldnames=fieldnames)
   #      amp_fixed_overview = list(reader)
   
   data_detailed = []
   for row in amp_overview:
      try:
         data_detailed.append([
            row['projectName'], 
            row['className'],
            row['numberOfTestMethodsBeforeShreding'], #before shreding
            row['targetLoc'], 
            row['mutationScoreBefore'], 
            row['n_notCoveredInOriginal'],
            row['n_mutantsInOriginal'], 
            row['testClassTimeToRunInMillis'], 
            workflow_params['mode'], # mode
            folder, # iteration-> I use folder name that includes the run Id
            row['status'], 
            row['numberOfOriginalTestMethods'], #after shreding
            row['timeBudgetFinished'], 
            row['n_amplifiedMethods'],
            row['n_newCovered'], 
            row['mutationScoreAfter'],
            row['mutationImprove'],
            row['numberOfProcessedMethods'], 
            row['duplicateMutants'],
            row['timeTotal']
         ])
      except KeyError as x:
         print(row)
         print(row)
         raise x
   header_row_out = ['projectName', 'className', 'numberOfTestMethodsBeforeShreding', 'targetLoc', 
            'mutationScoreBefore', 
            'n_notCoveredInOriginal',
            'n_mutantsInOriginal', 
            'testClassTimeToRunInMillis', 
            'mode', 
            'folder', 
            'status', 
            'numberOfOriginalTestMethods',
            'timeBudgetFinished', 
            'n_amplifiedMethods',
            'n_newCovered', 
            'mutationScoreAfter',
            'mutationImprove',
            'numberOfProcessedMethods',
            'duplicateMutants',
            'timeTotal']
   print(','.join(header_row_out))
   for row in data_detailed:
      print(','.join([str(x) for x in row]))

def createZipFilesNames(basename):
   # smallAmp-logs-Roassal3-run66
   nameParts = basename.split('.')[0].split('-')
   if len(nameParts)<4:
      sys.exit('I expect a file name like: smallAmp-logs-Roassal3-run66.zip')
   names = ['-'.join([nameParts[0], str, nameParts[2], nameParts[3]]) + '.zip' for str in ['logs', 'overview', 'results']]
   return names

def createFolderName(zipFiles):
   nameParts = zipFiles[0].split('.')[0].split('-')
   return '-'.join([ nameParts[2], nameParts[3] ])

def extractZips(_aZipAddress):
   dirname = os.path.dirname(_aZipAddress)
   basename = os.path.basename(_aZipAddress)
   zipFiles = createZipFilesNames(basename)
   zipFiles = [name for name in zipFiles if os.path.exists(os.path.join(dirname,name))]
   print('Available zip files: ', zipFiles)
   if len(zipFiles) < 3:
      sys.exit('There should be at least 3 zip files (overview, logs, results)')
   folder = os.path.join(dirname, createFolderName(zipFiles) )
   if os.path.exists(folder):
      sys.exit('Folder exist, run script using the folder name or delete it manually. Folder=' + folder)
   print('creating folder= ' + folder)
   os.mkdir(folder)
   
   for basename_zip in zipFiles:
      with zipfile.ZipFile(os.path.join(dirname, basename_zip), 'r') as zip_ref:
         zip_ref.extractall(folder)
      os.remove(os.path.join(dirname, basename_zip))
   
   print('zips are extracted: '+ str(os.listdir(folder)) )
   for zipName in os.listdir(folder):
      if zipName.endswith('.zip'):
         with zipfile.ZipFile(os.path.join(folder, zipName), 'r') as zip_ref:
            zip_ref.extractall(folder)
         os.remove(os.path.join(folder, zipName))
   print('internal zips are also extracted: '+ str(os.listdir(folder)) )
   return folder

def main():
   if len(sys.argv) < 2:
      sys.exit("I need a folder name")
   _arg = sys.argv[1]
   if not _arg.startswith(tmpDir):
      sys.exit("Move zip files in tmp dir.")
   if _arg.endswith('.zip'):
      print('extracting the files.')
      folder = extractZips(_arg)
      print('Afterwards, use the folder name to refer to this result. folder name=' + folder)
   else:
      folder = _arg
   process_folder(folder)

if __name__ == "__main__":
   print('run me from interactive python3')
   # main()
