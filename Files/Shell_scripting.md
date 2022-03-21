# Introduction to Shell Scripting

In the last tutorial we introduced the Bourne shell, a command line interpreter which allows us to execute commands, manipulate text and files and execute programs within a UNIX/LINUX system. However the shell also has excellent scripting features which allows us to automate tasks that would otherwise require us to enter multiple commands in succession and build and save analytical pipelines for use in multiple projects and tasks.

In its purest form a shell script is a list of commands listed in the order in which they are executed for executing a specific task or analysis. Shell scripts can be used to manipulate files, execute different programs in succession and print outputs such as texts, tables, figures etc. The shell is also a comprehensive programming language which we can use to write novel functions and programs but this feature of shell scripting is out of the scope of this tutorial. Here we will mostly concentrate on the ability of shell scripts to build and save analytical pipelines.

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
