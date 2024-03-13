# Introduction to Shell Scripting

In the last tutorial we introduced the Bourne shell, a command line interpreter which allows us to execute commands, manipulate text and files and execute programs within a UNIX/LINUX system. However the shell also has excellent scripting features which allows us to automate tasks that would otherwise require us to enter multiple commands in succession and build and save analytical pipelines for use in multiple projects and tasks.

In its purest form a shell script is a list of commands listed in the order in which they are executed for executing a specific task or analysis. Shell scripts can be used to manipulate files, execute different programs in succession and print outputs such as texts, tables, figures etc. The shell is also a comprehensive programming language which we can use to write novel functions and programs but this feature of shell scripting is out of the scope of this tutorial. Here we will mostly concentrate on the ability of shell scripts to build and save analytical pipelines.


- [A Basic Script](https://github.com/iksaglam/MBGE_411/blob/main/Files/Shell_scripting.md#a-basic-script)
- [Using Variables](https://github.com/iksaglam/MBGE_411/blob/main/Files/Shell_scripting.md#using-variables)
  - [Calling Values](https://github.com/iksaglam/MBGE_411/blob/main/Files/Shell_scripting.md#calling-values)
  - [Special Variables in UNIX/LINUX](https://github.com/iksaglam/MBGE_411/blob/main/Files/Shell_scripting.md#special-variables-n-in-unixlinux)
- [Example Scripts](https://github.com/iksaglam/MBGE_411/blob/main/Files/Shell_scripting.md#example-scripts)
  - [Script version 1](https://github.com/iksaglam/MBGE_411/blob/main/Files/Shell_scripting.md#script-version-1)
  - [Script version 2](https://github.com/iksaglam/MBGE_411/blob/main/Files/Shell_scripting.md#script-version-2)
  - [Script version 3](https://github.com/iksaglam/MBGE_411/blob/main/Files/Shell_scripting.md#script-version-3)


## A Basic Script

All shell scripts end with the extension `.sh` and start with a `shebang` construct in the first line. The `shebang` construct is an essential line which informs the machine which `language` should be used to interpret the script.

We will now create an example shell script called `welcome.sh`. To write a shell script we need to use a text editor. In this case we will be using `vim` as our editor of choice but feel free to use any text editor you wish. 

We first create an empty text file called `welcome.sh` using the following command

```Bash
vim welcome.sh
```

Next to the empty text we will add the `shebang` line. Don’t forget to enter `insert mode` in `vim` by hitting the `i` key.

```Bash
#!/bin/bash 
```

This is called a `shebang` line because the `#` symbol is called a `hash`, and the `!` symbol is called a `bang`. This line simply informs our machine to use the `bash` language located within the `bin` directory to interpret the script. For example, if the language of the script was written in `python` we would instead start our script with a `shebang` line designating the python environment.

```Bash
#!/usr/bin/env python
```

Next we will want to populate our script with the commands we want to execute. Each new command should be written in a different line in the order in which we want them executed. For example, let us write a simple script which will output a custom welcome message to the screen depending on the input we give it.

```Bash
#!/bin/bash 
echo "What is your name?"
read PERSON
echo "Welcome back $PERSON :)"
```

The above script passes the message `what is your name?` on to the screen and then uses the `read` command to assign the answer we type into a variable called `PERSON` and then prints the value of this variable on to our screen. Here the `$` sign stands for `get`.

We can execute the script using the following command:

```Bash
sh welcome.sh
```

This is what the output would look like:

```Bash
sh welcome.sh 
What is your name?
Ismail
Welcome back Ismail :)
```

If we want we can also make our script executable:

```Bash
chmod +x welcome.sh
```

This way we can execute the script directly without having to designate what type of script it is.

```Bash
./welcome.sh
What is your name?
Ismail
Welcome back Ismail :)
```

## Using Variables

In the above example as you might have noticed we passed a `variable` directly from the screen to our `shell` script. A `variable` is a character string to which a `value` is assigned. This value can be a specific `text`, `number`, `filename`, `program`, or any type of input you can think of.

We can define a `variable` within the `shell` like below:

```Bash
my_variable_name=variable_value
```

Some examples:
```Bash
LOC=“Sariyer”
infile=uvarovi.fa
var1=300
```

#### Calling Values

To call a certain value, text or file stored as a variable, all we have to do is `prefix` the name of the variable with the dollar sign `$`.

The following script `loc.sh` will call for the variable defined as `LOC` and `print` it onto the screen (i.e. `STDOUT`)

```Bash
#!/bin/bash
LOC=“Sariyer Istanbul”
echo $LOC
```
```Bash
sh loc.sh
Sariyer Istanbul
```

#### Special Variables `$n` in UNIX/LINUX

These variables allow us to pass `argument`s on to our scripts based on the order in which they are invoked. In this case `n` is a positive decimal number corresponding to the `position of the argument(s)` that directly come after our script in the command prompt (the first argument is `$1`, the second argument is `$2`, etc.).

For example the following script `personal.loc.sh` uses command line arguments `NAME` and `LOC` to write the location of a person on to the screen

```Bash
#!/bin/bash
NAME=$1
LOC=$2
echo $NAME is located in $LOC 
```
```Bash
sh personal.loc.sh  Ismail Sariyer
Ismail is located in Sariyer
```
## Example Scripts

#### Script version 1

Now let us write a `shell script` for a specific task. In this case we have several reference genomes within the directory `~/hpc_run/2024SpringMBGE411/week02_tutorial/genomes`. We want to `index` all of these genomes using the program `samtools` and from the indexed file count the number of chromosomes that are over 1 million base pairs in length.

To execute this procedure for a single genome we could write a script like below:

```Bash
#!/bin/bash
infile=$1
filter=$2
samtools faidx ${infile}.fa
awk "\$2 > $filter" ${infile}.fa.fai | wc -l | sed "s/^/$infile\t/" - > count.txt
```

Here the script requires `two arguments`: the `name of the genome` and the `number of basepairs` to filter. We can execute this script for a given genome (for example `Gryllus bimaculatus`) as below:

```Bash
sh index.ref.filter_v1.sh Gryllus.bimaculatus 1000000
```

The output is an `index` file ending in `.fai` (`Gryllus.bimaculatus.fa.fai`) and a file called `count.txt` which contains the number of chromosomes over 1000000 in Gryllus bimaculatus

#### Script version 2

Although the abobe is a useful script it requires us to execute the script for each genome separately. However, it would be much more practical if our script could `loop` through a `list` of genomes and execute the same command for each automatically. This way we could do the same job in only a single command.

To achieve this we need to slightly modify our script by adding in a `loop function` and also changing our first argument from a `single item` into a `list`. Below is a script that does exactly this:

```Bash
#!/bin/bash

infile=$1 # file containing list of individuals
filter=$2 # value to filter

for species in `cat $infile`;
do
samtools faidx ${species}.fa
awk "\$2 > $filter" ${species}.fa.fai | wc -l | sed "s/^/$species\t/" - >> count_indexes.txt
done
```

Before executing the script we need to prepare a list containing the names of all of our genomes. We can do this with a simple command as below:

```Bash
ls *.fa | cut -d'.' -f1-2 > species.list
```

Now that our list is ready we can execute the script like below:

```Bash
sh index.ref.filter_v2.sh species.list 1000000
```

#### Script version 3

The above is a pretty good script. The only thing missing is that we are not really using the available resources all that efficiently. One good thing about having access to a supercluster is that we can use dedicated computational nodes for running our analysis. This not only allows us to use much more processing power and memory, it also helps us to run multiple jobs at the same time or work on different stuff while our analysis is running.

To `submit` our job to a dedicated computing cluster we have to again slightly modify our script. Specifically we will have to designate in our script the resources we wish to use. We can achieve this by using the `$SBATCH` command:

```Bash
#!/bin/bash

#SBATCH --job-name=count_index
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --partition=short
#SBATCH --time=30

infile=$1 # file containing list of individuals
filter=$2 # value to filter

for species in `cat $infile`;
do
samtools faidx ${species}.fa
awk "\$2 > $filter" ${species}.fa.fai | wc -l | sed "s/^/$species\t/" - >> count_indexes.txt
done
```

Here we have requested a single core `#SBATCH --ntasks=1` from a single node `#SBATCH --nodes=1` and 16G of memory `#SBATCH --mem=16G`. We have also designated that we wish to use the short partition `#SBATCH --partition=short` and have requested 30 minutes of wall time `#SBATCH --time=30`. We have also named our job `#SBATCH --job-name=count_index` because why not! 

Now all we have to do is submit our script to the que where it will automatically get placed in an appropriate computing node. To do this we use the `sbatch` command.

```Bash
sbatch index.ref.filter_v3.sh species.list 1000000 
```

We can check the details of our job and fing information about how long it has been running, which specific node and partition it is on using the squeue command.

```Bash
squeue -u mbge411
             JOBID  PARTITION   NAME      USER      ST      TIME      NODES   NODELIST(REASON)
           2811657  short       count_in  mbge411   R       0:02      1       it03
```

This command will also produce a slurm file `slurm-${JOBID}.out` which is a log file of all activities going on and is very useful for parsing errors. Always make sure to check the output in your slurm file before moving on to the next or before looking at your results!


