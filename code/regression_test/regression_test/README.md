# Regression Test for liblinear

### Requirements

git

Optional: gnuplot, matlab

### Preparation

First

```shell
$ ./get.sh [OLD_HASH]
```
default OLD_HASH is HEAD^

Then 'old/' and 'new/' will be created:
old is the dir contain the version of OLD_HASH
new is the dir contain the version of HEAD

### Execution

Then running

```shell
$ ./regression.py old new [cmds] [results/]
```

will perform regression tests, and the terminal outputs will also be stored in `results/log.txt`. Detailed results will be put under `results/output/`.

### Results

The results will present in the following format
testID status cmds

testID: represent the index of the commands.
status: there are three type of status. 
	1. pass: all result is the same. 
	2. diff: give the file that changed (you could check the 'results/output/' in detail).
	3. skip: the command cannot test in this enviroment.
	4. error: the command exit abnormally.
cmds: the test command.
