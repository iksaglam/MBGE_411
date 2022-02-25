
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

### File Manipulation: Making a directory

Creating a directory is pretty easy. The command we use is `mkdir` which stands for Make Directory. Using this command let us create a new directory under our own name

```Bash
mbge411@login02:~$ ls 
bin  scripts  Spring2019  Spring2020  week01_tutorial
mbge411@login02:~$ mkdir ismail
mbge411@login02:~$ ls
bin  ismail  scripts  Spring2019  Spring2020  week01_tutorial
```


### File Manipulation: Removing a Directory

Just like creating a directory, removing a directory is pretty straightforward and depends on a single command. In this case the command is `rmdir`, which stands for remove directory. However, be very careful when using the `rmdir` command because there is no undo when it comes to the command line. So if you remove or delete a directory it will be gone forever. For example, let us first create a new directory called `temp` under `ismail/` and then remove it.

```Bash
mbge411@login02:~$ mkdir ~/ismail/temp
mbge411@login02:~$ cd ismail
mbge411@login02:~/ismail$ ls
temp
mbge411@login02:~/ismail$ rmdir temp
mbge411@login02:~/ismail$ ls

```

### File Manipulation: Creating a Blank File

A lot of commands that involve manipulating data in unix/linux have the nice feature that they will create a file automatically if we call upon a file that does not exist. We can use this to our advantage and create a blank file using the command `touch`. 

```Bash
mbge411@login02:~/ismail$ touch example1.txt
mbge411@login02:~/ismail$ touch example2.txt
mbge411@login02:~/ismail$ ls
example1.txt example2.txt
```

### File Manipulation: Copying a File or Directory

There are many instances in which we might want to create a duplicate of a certain file or copy an existing file into another directory. The command for this is `c`p which stands for copy. Using `cp` is straightforward but note that copying directories requires the use of the `-r` option. 

```Bash
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
```

### File Manipulation: Moving a File or Directory

To move a file we use the command `mv`, which is short for move. It is very similar to `cp` but can be directly used on directories and does not require the `-r` option.

```Bash
mbge411@login02:~/ismail$ ls
example1.txt example2.txt example/ example2/
mbge411@login02:~/ismail$ mv example/ example3/
mbge411@login02:~/ismail$ ls
example1.txt example2.txt example2/ example3/
mbge411@login02:~/ismail$ mv example2/ example3/example4
mbge411@login02:~/ismail$ cd example3/
mbge411@login02:~/ismail/example3$ ls
example4/
```

We can also use the `mv` command to also rename files or directories. For example:

```Bash
mbge411@login02:~/ismail$ ls
example1.txt example2.txt example3/
mbge411@login02:~/ismail$ mv example1.txt example3.txt
mbge411@login02:~/ismail$ ls
example2.txt example3.txt example3/
```

### File Manipulation: Removing a file

The command to remove or delete a file is `rm` which stands for remove. As with `rmdir` the process is non-reversible so be very careful before using it.

```Bash
mbge411@login02:~/ismail$ ls
example2.txt example3.txt example3/
mbge411@login02:~/ismail$ ls
rm example2.txt
mbge411@login02:~/ismail$ ls
example3.txt example3/
mbge411@login02:~/ismail$ rm -r example3/
mbge411@login02:~/ismail$ ls
example3.txt
```

Note that this time to remove the directory `example3/` we used the command rm `-r` and not `rmdir`. This is because `rmdir` only works with empty directories. If directories have files or subdirectories within them `rmdir` will not work.


### Wildcards

Above we learned some useful commands to manipulate files and directories. However, as you might have noticed we manipulated files or directories one by one. When we have a lot of files this might not be practical and having a way to manipulate several files at once would be really useful.

This is where wildcards come in. Wildcards are a set of symbols that allow you to create a pattern defining a set of files or directories. There is not enough space here to go over all wildcards but some basic wildcards are given below:

- `*` -> means match all characters
- `?` -> means match any single character
- `[ ]` -> represents a range of characters

Let us start with the wildcard you will use most, `*`. In the examples below we will `list` everything beginning with `e`.

```Bash
mbge411@login02:~$ cd week01_tutorial/wildcards
mbge411@login02:~/week01_tutorial/wildcards$ ls
example1.txt example2.txt example3.txt example4.txt example5.txt foo1.jpg foo2.jpg foo3.jpg foo4.jpg gmail1.thml gmail2.html gmail3.html 
mbge411@login02:~/week01_tutorial/wildcards$ ls e*
example1.txt example2.txt example3.txt example4.txt example5.txt
```

We can use `*` in a number of other ways. For example, let us use this wildcard to list only image `jpg` files.

```Bash
mbge411@login02:~/week01_tutorial/wildcards$ ls *.jpg
foo1.jpg foo2.jpg foo3.jpg foo4.jpg
```

Or to list all files containing the letter `m`:

```Bash
mbge411@login02:~/week01_tutorial/wildcards$ ls *m*
example1.txt example2.txt example3.txt example4.txt example5.txt gmail1.thml gmail2.html gmail3.html
```

Now let us look at how we can use the wildcard `?`. In the example below we wish to list all `foo` files.

```Bash
mbge411@login02:~/week01_tutorial/wildcards$ ls foo?.jpg
foo1.jpg foo2.jpg foo3.jpg foo4.jpg
```

We can also use wildcards together to come up with even more detailed filters. For example, let us list all files that have three letter extensions.

```Bash
mbge411@login02:~/week01_tutorial/wildcards$ ls *.???
example1.txt example2.txt example3.txt example4.txt example5.txt foo1.jpg foo2.jpg foo3.jpg foo4.jpg
```

Or let us list only those files whose third letter is “a”

```Bash
mbge411@login02:~/week01_tutorial/wildcards$ ls ??a*
example1.txt example2.txt example3.txt example4.txt example5.txt gmail1.thml gmail2.html gmail3.html
```

Finally let us introduce the range `[ ]` operator. This operator allows you to limit your calls to a certain subset of characters. For example let us list all files that either begin with an `e` or `g`

```Bash
mbge411@login02:~/week01_tutorial/wildcards$ ls [eg]*
example1.txt example2.txt example3.txt example4.txt example5.txt gmail1.thml gmail2.html gmail3.html
```

Or let us list all files that contain numbers between `1` and `3`:

```Bash
mbge411@login02:~/week01_tutorial/wildcards$ ls *[1-3]*
example1.txt example2.txt example3.txt foo1.jpg foo2.jpg foo3.jpg gmail1.thml gmail2.html gmail3.html
```

### Viewing and editing files

In the previous sections we created and manipulated files that were blank. Of course in reality we will not be working with blank files but with read data (i.e. files that contain information). An important part of our workflow will entail working with files containing information and creating and writing new information to new files. 

A good way to write new information or content into blank files is to use a command line text editor (similar to text editors like `Notepad` or `TextWrangler` unique to Windows and Mac respectively).

One such command line text editor is `Vim`. There certainly are other text editors you might wish to try (like `Emacs` or `Nano`) but for now we will stick with `Vim` because it’s the one I like and hence the best!!! The most important aspect of working with `Vim` (or any other command line text editor) is that everything is done through the `terminal` so you will have to say goodbye to the `mouse` and start using the `up/down` and `left/right` arrow keys to move!!!

There are two modes in `Vim`. `Insert` mode and `Edit` mode. In `insert` mode you may input or enter content into the file (in other words you can write to the file). In `edit` mode you can move around the file, perform actions such as delete, copy, search, replace, saving etc. 

To open up any file with `Vim` we only need the `vim` command followed by the `file` we would like to open. Like all linux tools if the file you specified does not exist, your computer will create one for you automatically.

```Bash
mbge411@login02:~$ vim myfile.txt
```

When you run this command it will open up the file. Since no such file existed, the file will be empty and what you see will look something like this. 

```Bash
~
~
~
~
"Myfile.txt" [New File]
```

`Vim` always starts in `edit` mode so the first thing you will want to do is switch to `insert` mode so we can write something to our text file. To enter `insert` mode all we have to do is hit the `i` key. Once we have done this we can check that we are indeed in `insert` mode by looking at the left hand corner, which should now display `--INSERT--`

```Bash
~
~
~
~
~
-- INSERT --
```

Now that we are in `insert` mode go ahead and type in something. For example I can type in my name

```Bash
ismail
~
~
~
~
-- INSERT --
```

When finished writing you can hit `Esc` (escape) to exit `insert` mode and go back to `edit` mode. Now that we have finished editing our text file, we would naturally want to `save` and `exit`. There are several ways we can accomplish this in `Vim` which are given below. All will get the job done but make sure you are in `edit` mode before trying to use them!

```Bash
ZZ (capitals!!) - Save and exit
:q! - discard all changes since the last save and exit
:w - save file but don't exit
:wq - save and exit
:x - save and exit
```

Of course `Vim` has a host of other commands for doing lots of other stuff. [Here](https://vim.rtorr.com/) is a good cheat sheet!!

### Other ways to view files

`Vim` will allow us to edit files, but sometimes we simply want to view the contents of a file without worrying about editing. In these instances it might be more convenient to use the `cat` or `less` commands.

The `cat` command stands for `concatenate` and it’s main purpose is to join two files together. However in its most simple form it can be used to display the contents of any file on to the screen. For example let us run the `cat` command with our newly created text file:

```Bash
mbge411@login02:~$ cat myfile.txt 
Ismail
```

The `cat` command is useful when we want to view a small file, but if the file we want to view is large with hundreds or thousands of lines then writing all of its content onto the screen is not practical. For such large files a better suited command is `less`.

```Bash
mbge411@login02:~$ less mylongfile.txt
```

`less` allows you to move up and down the file using `arrow keys`. You can even jump a whole page using the `SpaceBar` or move back a page by pressing `b`. You can jump to the end of the file by pressing `Shift+g` or go back to the top of the file by pressing `gg`. 

### Filters

The `Unix/Linux` command line is such a powerful tool for filtering data that once you have a basic understanding of it, it will be hard to go back to GUI spreadsheets (i.e. `Excel` and the like). Filters in `Unix/Linux` are specific commands that accept textual data and transform it in a way more suited to our needs or what we wish to see.

We will introduce the usage of these commands by direct example and manipulation of the file `isophya.csv`. Also, please remember that almost all commands introduced below have multiple options so make sure to check out their `man` pages.

Before moving on to specific commands let us first view our example file using the `cat` command.

```Bash
mbge411@login02:~$ cd week01_tutorial/filters/
mbge411@login02:~/week01_tutorial/filters$ ls
isophya.csv uarctos.metadata
mbge411@login02:~/week01_tutorial/filters$ cat isophya.csv
IST06	Dark	450	-0.94	65.9
IST06	Dark	450	-0.94	65.9
PLVYL	Dark	850	-0.77	73.7
PLVYL	Dark	850	-0.77	73.7
IST13	Dark	1000	-0.86	75.66
IST13	Dark	1000	-0.86	75.66
PIKNK	Pale	1300	0.58	78
VRCNK	Pale	2000	0.81	76.92
KALEK	Pale	2100	0.05	76.74
KALEK	Pale	2100	0.05	76.74
KALEK	Pale	2100	0.05	76.74
CKLYU	Pale	2300	0.00	75.81
CKLYU	Pale	2300	0.00	75.81
CKLYU	Pale	2300	0.00	75.81
```

