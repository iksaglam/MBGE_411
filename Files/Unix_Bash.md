
# Introduction to Linux and Bash 

### Opening a Terminal

Opening a terminal is easy but will differ from system to system. Below are a few places to start looking depending on which system you are on.

- On a Linux machine you will find it in `Applications -> System or Applications -> Utilities`.
- On a Mac you will find it under `Applications -> Utilities`.
- On a Windows machine the terminal is called a command prompt and you can find this under `Start -> Program Files -> Accessories -> Command Prompt`. However, Windows does not have a native bash shell so in Windows 10 you will need to separately install [Bash](https://www.windowscentral.com/how-install-bash-shell-command-line-windows-10) or if you are using an older version of windows you will need to install an external SSH client. [Putty](http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html) or [Cygwin](https://cygwin.com/install.html) are rather good ones.

### Connect to a server using SSH

SSH is a protocol through which you can access your remote server and run shell commands. SSH is encrypted with Secure Sockets Layer (SSL), which makes it difficult for these communications to be intercepted and read.

Using the Internet Protocol (IP) address and password you can connect (or log in) to your remote server following the below command:

```Bash
ssh mbge411@login.kuacc.ku.edu.tr
```


The system will prompt you to enter a password for the account you are trying to connect to. Type in the below password and hit enter.

```Bash
Ap9r3LqVVm
```

Notice that as you type, none of the characters will be displayed on the screen. This is normal so be very careful to type correctly.


### The Shell, Bash

Once logged in (or connected) to our remote server we will be greeted by the shell interface. This is a part of the operating system that defines how you will interact with your computer and how it will behave and look after running (or executing) commands. There are various shells available but the most common and the one linux computers (such as our remote server) uses is called bash which stands for Bourne again shell. 

Let us test this and see if indeed we are using bash as our shell. For this we will use a command called echo which is used to display messages directly on our screen.

```Bash
mbge411@login02:~$ echo $SHELL
/bin/bash
```

As we can see this command prints out something ending in bash and verifies that we are indeed using the bash shell.

### Navigating the Shell: Where are we?

Once we have logged on to our server the first thing we want to know is where we are. To do this we will enter the command pwd which stands for Print Working Directory.

```Bash
mbge411@login02:~$ pwd
/kuacc/users/mbge411
````

A lot of the commands we are executing will depend on where we are or where our target files are. Therefore, keeping track of where we are using the pwd command is good practice.

### Navigating the Shell: What is in our directory? 

Once we know where we are the next thing we will want to do is list the contents of the directory we are in. The command to do this is `ls` which as you might have guessed stands for list.

```Bash
mbge411@login02:~$ ls
bin  course_content  scripts  Spring2021
```

Whereas pwd just ran by itself, ls is a bit more powerful and can run with arguments which enable us to see more details concerning our files. For example lets run `ls` with the argument `-l`.

```Bash
mbge411@login02:~$ ls -l
total 2
drwxr-s--- 10 mbge411 gde  9 May 14  2019 bin
drwxr-s---  2 mbge411 gde  0 Feb 12  2019 scripts
drwxr-sr-x 14 mbge411 gde 12 Jan 31 16:40 Spring2019
drwxr-sr-x  2 mbge411 gde  0 Jan 31 16:55 Spring2020
drwxr-sr-x  4 mbge411 gde  2 Jan 31 16:37 week01_tutorial
```


`ls -l` stands for long listing which as you can see gives us a lot more information about our files. Let’s break them down.

- The first character indicates whether it is a normal file (-) or directory (-d)
- Next 9 characters are permissions to our file or directory (we will not worry about them now).
- The next field is the number of blocks (again nothing that concerns us).
- The next field is the owner of the file or directory (mbge411 in this case).
- The next field is the group the file or directory belongs to (gde in this case).
- The next field is the file size (0, since they are empty).
- The next field is the file modification time.
- Lastly we have the actual name of the file or directory.

The `ls` command has many more arguments for a variety of purposes. We don’t have time to go over all of them here. To view all possible usage of `ls` we can open up it’s manual by typing the command `ls man`, where man stands for manual. This is a good way of learning about any other command in the shell interface. Therefore, whenever we are using a new command it is always a good idea to hit man and learn all that command has to offer. Now type ls man in your terminal, hit enter and see what happens.

### Navigating the Shell: Paths 

Whenever we refer to a file or directory on the command line, we are in reality referring to a path. A path is a set of directions to get to a particular file or directory on the system.
There are two types of paths, `absolute` and `relative`. We can use whichever type of path we prefer as the system will direct us to the same location.

The file system under linux (or Unix) is a hierarchical structure. At the very top is the root directory, denoted by a single slash `-/-`. Below the root directory are many subdirectories. These subdirectories have subdirectories of their own etc. Files may reside in any one of these directories.

`Absolute paths` specify a location (file or directory) in relation to the root directory and always begin with a forward slash `-/-`.

`Relative paths` specify a location (file or directory) in relation to where we currently are in the system. They will not begin with a slash.

Paths also can be designated using shorthand notations.
- `~ (tilde)`**:** This is a shorthand notation for your home directory. For example if our home directory is `/kuacc/users/mbge411`, then we could refer to the directory Sprin2021 with the path  `/kuacc/users/mbge411/Spring2021` or simply as `~/Spring2021`.
- `. (dot)`**:** This is a shorthand notation for the current directory. In the above example we referred to this directory using it’s full path. It could also be written as `./Spring2021`.
- `.. (dotdot)`**:** This is a shorthand notation to the parent directory or to the directory one step above of the current directory we are in. For example, if we were in the path `/kuacc/users/mbge411/Spring2021` and we typed in `ls ../` this would give us a listing of our home directory `/kuacc/users/mbge411/`

It should be obvious by now that we can refer to a specific location in multiple ways. The important thing to note here is that there is no one correct way. Each way is equally valid and whichever you choose to use mostly will depend on how you like to do things. Below are a few examples:

```Bash
mbge411@login02:~$ pwd
/kuacc/users/mbge411
mbge411@login02:~$ ls ~/week01_tutorial
filters  wildcards
mbge411@login02:~$ ls ./week01_tutorial
filters  wildcards
mbge411@login02:~$ ls /kuacc/users/mbge411/week01_tutorial
filters  wildcards
mbge411@login02:~$ ls ../../
admin  ALU_LEFT  bashconfig  images  kuacc  users
mbge411@login02:~$ ls /
ansys_inc  bin  boot  build  dev  etc  ext  home  home2  kuacc ...
```

Try playing around with the above commands to get a feel for how paths work and what they mean.

### Navigating the Shell: Moving around 

Now that we know what directories and paths are, the next thing we would want to do is move between them. To move around the system we use a simple command called `cd` which as you might have guessed stands for change directory. To use it you simply type `cd` and the path to where you want to go: `cd <location>`

If you do not specify a path and just type `cd` on your terminal, no matter where you are, it will take you back to your `home directory`. To change directories with the `cd` command we can supply it with  an `absolute path` or `relative path`. Either will work and get us to where we want to go. Below are a few examples.

```Bash
mbge411@login02:~$ pwd
/kuacc/users/mbge411
mbge411@login02:~$ ls 
bin  scripts  Spring2019  Spring2020  week01_tutorial
mbge411@login02:~$ cd week01_tutorial
mbge411@login02:~/week01_tutorial$ ls
filters  wildcards
mbge411@login02:~/week01_tutorial$ cd ~/bin
mbge411@login02:~/bin$ ls
FastQC  fastx-toolkit  gffcompare-0.10.6.Linux_x86_64 ...
mbge411@login02:~/bin$ cd ../../../
mbge411@login02:/kuacc$ ls
apps  etc  jobscripts  kadmin  users
mbge411@login02:/kuacc$ cd
mbge411@login02:~$ ls
bin  scripts  Spring2019  Spring2020  week01_tutorial
```

File Manipulation: Making a directory

Creating a directory is pretty easy. The command we use is mkdir which stands for Make Directory. Using this command let us create a new directory under our own name

mbge411@login02:~$ ls 
bin  scripts  Spring2019  Spring2020  week01_tutorial
mbge411@login02:~$ mkdir ismail
mbge411@login02:~$ ls
bin  ismail  scripts  Spring2019  Spring2020  week01_tutorial
mbge411@login02:~$



File Manipulation: Removing a Directory

Just like creating a directory, removing a directory is pretty straightforward and depends on a single command. In this case the command is rmdir, which stands for remove directory. However, be very careful when using the rmdir command because there is no undo when it comes to the command line. So if you remove or delete a directory it will be gone forever. For example, let us first create a new directory called temp under ismail and then remove it.

mbge411@login02:~$ mkdir ~/ismail/temp
mbge411@login02:~$ cd ismail
mbge411@login02:~/ismail$ ls
temp
mbge411@login02:~/ismail$ rmdir temp
mbge411@login02:~/ismail$ ls

mbge411@login02:~/ismail$




File Manipulation: Creating a Blank File

A lot of commands that involve manipulating data in linux have the nice feature that they will create a file automatically if we call upon a file that does not exist. We can use this to our advantage and create a blank file using the command touch. 

mbge411@login02:~/ismail$ touch example1.txt
mbge411@login02:~/ismail$ touch example2.txt
mbge411@login02:~/ismail$ ls
example1.txt example2.txt
mbge411@login02:~/ismail$



File Manipulation: Copying a File or Directory

There are many instances in which we might want to create a duplicate of a certain file or copy an existing file into another directory. The command for this is cp which stands for copy. Using cp is straightforward but note that copying directories requires the use of the -r option. 

mbge411@login02:~/ismail$ ls
example1.txt example2.txt
mbge411@login02:~/ismail$ cp example1.txt example1_copy.txt
mbge411@login02:~/ismail$ ls
example1_copy.txt example1.txt example2.txt
mbge411@login02:~/ismail$ cp example1.txt ../
mbge411@login02:~/ismail$ ls ../
bin  ismail  example1.txt  scripts  Spring2019  Spring2020 ...
mbge411@login02:~/ismail$ mkdir example
mbge411@login02:~/ismail$ ls
example1.txt example2.txt example/
mbge411@login02:~/ismail$ cp -r example example2
mbge411@login02:~/ismail$ ls
example1.txt example2.txt example/ example2/



File Manipulation: Moving a File or Directory

To move a file we use the command mv, which is short for move. It is very similar to cp but can be directly used on directories and does not require the -r option.

mbge411@login02:~/ismail$ ls
example1.txt example2.txt example/ example2/
mbge411@login02:~/ismail$ mv example/ example3/
mbge411@login02:~/ismail$ ls
example1.txt example2.txt example2/ example3/
mbge411@login02:~/ismail$ mv example2/ example3/example4
mbge411@login02:~/ismail$ cd example3/
mbge411@login02:~/ismail/example3$ ls
example4/
mbge411@login02:~/ismail/example3


We can also use the mv command to also rename files or directories. For example:

mbge411@login02:~/ismail$ ls
example1.txt example2.txt example3/
mbge411@login02:~/ismail$ mv example1.txt example3.txt
mbge411@login02:~/ismail$ ls
example2.txt example3.txt example3/



File Manipulation: Removing a file

The command to remove or delete a file is rm which stands for remove. As with rmdir the process is non-reversible so be very careful before using it.

mbge411@login02:~/ismail$ ls
example2.txt example3.txt example3/
mbge411@login02:~/ismail$ ls
rm example2.txt
mbge411@login02:~/ismail$ ls
example3.txt example3/
mbge411@login02:~/ismail$ rm -r example3/
mbge411@login02:~/ismail$ ls
example3.txt


Note that this time to remove the directory example3/ we used the command rm -r and not rmdir. This is because rmdir only works with empty directories. If directories have files or subdirectories within them rmdir will not work.


