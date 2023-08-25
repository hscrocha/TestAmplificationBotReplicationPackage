# Test Amplification Bot for Pharo/Smalltalk (Reproduction Package)

Reproduction Package for our paper A Test Amplification Bot for PHARO/SMALLTALK (working submission)

## Abstract
Test amplification is exploiting the knowledge embedded in the existing test suite to strengthen it. A typical test amplification technique is transforming the initial tests to find new test methods that increase the mutation coverage.
Although past researches prove the benefits of this technique, they pay less attention to how to bring these tools to the everyday life of developers. 
Test amplifiers are generally complicated, and their execution takes a considerable time. 
The current state-of-the-art test amplifiers lack a mechanism to prioritize their task within a limited budget.
This paper integrates Small-Amp to Github-Actions to introduce a zero-touch test amplification.
In addition, we introduce a heuristic to prioritize test methods covering more alive mutants to increase the performance in the limited time.
Since Small-Amp amplifies the tests in a live system, Pharo, we also introduce a sandboxing mechanism to make it crash resilient.
We evaluated this approach by installing it on five open-source Pharo projects.
The experiment result shows that test amplification can be employed at a project level by integrating it into the build system.
It also shows that our prioritizing mechanism also successfully increases the tool performance from 6.28 mutant per class to 8.42 when execution is run out of time.
In conclusion, by integrating the mutation-testing-based test amplification tools in the continuous integration platforms, these tools are not a new burden on the shoulders of developers anymore.

