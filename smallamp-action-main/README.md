# SmallAmp Action

Use this action to run [SmallAmp](https://github.com/mabdi/small-amp/) on the Github Actions framework.

# How to use me

Look at our toy example project [SmallBank](https://github.com/mabdi/smalltalk-SmallBank).
You can Copy the workflow from [SmallAmpCI.yml](https://github.com/mabdi/smalltalk-SmallBank/blob/master/.github/workflows/SmallAmpCI.yml) and paste it in your project.

What you need to change in the workflow:

- **Jobs number**: 
To set the number of parallel jobs, you need to change the lengh of the array in `jobs.TestAmplification.strategy.matrix.portion`.
The content of array is not important. Just the length is important. 
Github actions supports up to 8 parallel jobs.

- **Metacello parameters**: 
Set `env.project_baseline`, `env.project_directory`, and `env.project_load`.
Setting `env.project_load` is optional, but other parameters are necessary.
Here is psudocode of how `Metacello` is called:
```
     b:= Metacello new
        baseline: {env.project_baseline};
        repository: 'tonel://', {env.clone_directory+ env.project_directory}.
     {env.project_load} ifNotNil: [ b load: {env.project_load} ] ifNil: [ b load ]
```
