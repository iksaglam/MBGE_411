
# Introduction to Unix and Bash 

- [Opening a Terminal](https://github.com/iksaglam/MBGE_411/blob/main/Files/Unix_Bash.md#opening-a-terminal)
- [Connect to a server using SSH](https://github.com/iksaglam/MBGE_411/blob/main/Files/Unix_Bash.md#connect-to-a-server-using-ssh)
- [The Shell: Bash](https://github.com/iksaglam/MBGE_411/blob/main/Files/Unix_Bash.md#the-shell-bash)
- [Navigating the Shell](https://github.com/iksaglam/MBGE_411/blob/main/Files/Unix_Bash.md#navigating-the-shell)
  - [Where are we?](https://github.com/iksaglam/MBGE_411/blob/main/Files/Unix_Bash.md#navigating-the-shell-where-are-we)
  - [What is in our directory?](https://github.com/iksaglam/MBGE_411/blob/main/Files/Unix_Bash.md#navigating-the-shell-what-is-in-our-directory)
  - [Paths](https://github.com/iksaglam/MBGE_411/blob/main/Files/Unix_Bash.md#navigating-the-shell-paths)
  - [Moving around](https://github.com/iksaglam/MBGE_411/blob/main/Files/Unix_Bash.md#navigating-the-shell-moving-around)
- [File Manipulation](https://github.com/iksaglam/MBGE_411/blob/main/Files/Unix_Bash.md#file-manipulation)
  - [Making a directory](https://github.com/iksaglam/MBGE_411/blob/main/Files/Unix_Bash.md#file-manipulation-making-a-directory)
  - [Removing a Directory](https://github.com/iksaglam/MBGE_411/blob/main/Files/Unix_Bash.md#file-manipulation-removing-a-directory)
  - [Creating a Blank File](https://github.com/iksaglam/MBGE_411/blob/main/Files/Unix_Bash.md#file-manipulation-creating-a-blank-file)
  - [Copying a File or Directory](https://github.com/iksaglam/MBGE_411/blob/main/Files/Unix_Bash.md#file-manipulation-copying-a-file-or-directory)
  - [Moving a File or Directory](https://github.com/iksaglam/MBGE_411/blob/main/Files/Unix_Bash.md#file-manipulation-moving-a-file-or-directory)
  - [Removing a file](https://github.com/iksaglam/MBGE_411/blob/main/Files/Unix_Bash.md#file-manipulation-removing-a-file)
- [Wildcards](https://github.com/iksaglam/MBGE_411/blob/main/Files/Unix_Bash.md#wildcards)
- [Viewing and editing files](https://github.com/iksaglam/MBGE_411/blob/main/Files/Unix_Bash.md#viewing-and-editing-files)
- [Filters](https://github.com/iksaglam/MBGE_411/blob/main/Files/Unix_Bash.md#filters)
  - [Head](https://github.com/iksaglam/MBGE_411/blob/main/Files/Unix_Bash.md#head)
  - [Tail](https://github.com/iksaglam/MBGE_411/blob/main/Files/Unix_Bash.md#tail)
  - [Sort](https://github.com/iksaglam/MBGE_411/blob/main/Files/Unix_Bash.md#sort)
  - [wc](https://github.com/iksaglam/MBGE_411/blob/main/Files/Unix_Bash.md#wc)
  - [cut](https://github.com/iksaglam/MBGE_411/blob/main/Files/Unix_Bash.md#cut)
  - [uniq](https://github.com/iksaglam/MBGE_411/blob/main/Files/Unix_Bash.md#uniq)
- [sed](https://github.com/iksaglam/MBGE_411/blob/main/Files/Unix_Bash.md#sed)
- [Grep & Regular Expressions](https://github.com/iksaglam/MBGE_411/blob/main/Files/Unix_Bash.md#grep--regular-expressions)
- [Piping & redirecting](https://github.com/iksaglam/MBGE_411/blob/main/Files/Unix_Bash.md#piping--redirecting)
- [awk](https://github.com/iksaglam/MBGE_411/blob/main/Files/Unix_Bash.md#awk)

#

## Opening a Terminal

Opening a terminal is easy but will differ from system to system. Below are a few places to start looking depending on which system you are on.

- On a Linux machine you will find it in `Applications -> System or Applications -> Utilities`.
- On a Mac you will find it under `Applications -> Utilities`.
- On a Windows machine the terminal is called a command prompt and you can find this under `Start -> Program Files -> Accessories -> Command Prompt`. However, Windows does not have a native bash shell so in Windows 10 and 11 you will have to install a Linux terminal separately. See instructions [here](https://webme.ie/how-to-install-ubuntu-terminal-for-windows-10/#:~:text=Ubuntu%20on%20Windows%20allows%20you,apt%20and%20more%20on%20Windows.&text=First%20you%20need%20to%20enable,panel%20and%20click%20on%20Programs.&text=Then%20check%20the%20box%20for%20Windows%20Subsystem%20for%20Linux%20and%20click%20ok.) or [here](https://www.windowscentral.com/how-install-bash-shell-command-line-windows-10). If you are using an older version of windows you will need to install an external SSH client. [Putty](http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html) or [Cygwin](https://cygwin.com/install.html) are rather good ones.

## Connect to a server using SSH

SSH is a protocol through which you can access your remote server and run shell commands. SSH is encrypted with Secure Sockets Layer (SSL), which makes it difficult for these communications to be intercepted and read.

Using the Internet Protocol (IP) address and password you can connect (or log in) to your remote server following the below command:

```Bash
ssh mbge411@login.kuacc.ku.edu.tr
```


The system will prompt you to enter a password for the account you are trying to connect to. Type in the below password and hit enter.

```Bash
GenomicsBio!
```

Notice that as you type, none of the characters will be displayed on the screen. This is normal so be very careful to type correctly.


## The Shell: Bash

Once logged in (or connected) to our remote server we will be greeted by the shell interface. This is a part of the operating system that defines how you will interact with your computer and how it will behave and look after running (or executing) commands. There are various shells available but the most common and the one linux computers (such as our remote server) uses is called bash which stands for Bourne again shell. 

Let us test this and see if indeed we are using bash as our shell. For this we will use a command called echo which is used to display messages directly on our screen.

```Bash
mbge411@login02:~$ echo $SHELL
/bin/bash
```

As we can see this command prints out something ending in bash and verifies that we are indeed using the bash shell.

## Navigating the Shell
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
## File Manipulation
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


## Wildcards

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

## Viewing and editing files

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

## Filters

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
### Head

`Head` is a command that prints out a certain number of lines of a given file. The default number of lines is 10, but we can modify this number using the option `-n`.

```Bash
mbge411@login02:~/week01_tutorial/filters$ head isophya.csv
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
```

```Bash
mbge411@login02:~/week01_tutorial/filters$ head -n 6 isophya.csv
IST06	Dark	450	-0.94	65.9
IST06	Dark	450	-0.94	65.9
PLVYL	Dark	850	-0.77	73.7
PLVYL	Dark	850	-0.77	73.7
IST13	Dark	1000	-0.86	75.66
IST13	Dark	1000	-0.86	75.66
```

### Tail

`Tail` is the opposite of head and prints out the last certain number of lines of a given file. Default number of lines is 10, but we can modify this number using the option `-n`.

```Bash
mbge411@login02:~/week01_tutorial/filters$ tail isophya.csv
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

```Bash
mbge411@login02:~/week01_tutorial/filters$ tail -n 3 isophya.csv
CKLYU	Pale	2300	0.00	75.81
CKLYU	Pale	2300	0.00	75.81
CKLYU	Pale	2300	0.00	75.81
```

You can also do reverse filtering using `tail`. For example, if we wanted to filter out the first 2 lines of our file we could use the following command:

```Bash
mbge411@login02:~/week01_tutorial/filters$ tail -n +3 isophya.csv
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

Here we are essentially telling our file to start printing from the top 3rd line.

### Sort

`Sort`, does exactly what it says and sorts files. The command by default will `sort alphabetically` and use the first column. But there are many options to modify the sorting process as we will see below.

```Bash
mbge411@login02:~/week01_tutorial/filters$ sort isophya.csv
CKLYU	Pale	2300	0.00	75.81
CKLYU	Pale	2300	0.00	75.81
CKLYU	Pale	2300	0.00	75.81
IST06	Dark	450	-0.94	65.9
IST06	Dark	450	-0.94	65.9
IST13	Dark	1000	-0.86	75.66
IST13	Dark	1000	-0.86	75.66
KALEK	Pale	2100	0.05	76.74
KALEK	Pale	2100	0.05	76.74
KALEK	Pale	2100	0.05	76.74
PIKNK	Pale	1300	0.58	78
PLVYL	Dark	850	-0.77	73.7
PLVYL	Dark	850	-0.77	73.7
VRCNK	Pale	2000	0.81	76.92
```

We can also `sort` using any column by designating the desired column with the `-k` option. For example if we want to sort by column 3:

```Bash
mbge411@login02:~/week01_tutorial/filters$ sort -k3 isophya.csv
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
IST06	Dark	450	-0.94	65.9
IST06	Dark	450	-0.94	65.9
PLVYL	Dark	850	-0.77	73.7
PLVYL	Dark	850	-0.77	73.7
```

Notice above that the `-k3` option sorted our file `alphabetically`. However, since column three is numeric `sorting numerically` would make more sense. We can do this by using the `-n` option.

```Bash
mbge411@login02:~/week01_tutorial/filters$ sort -nk3 isophya.csv
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

Lastly, perhaps we want to `sort` in the `reverse order` (i.e. from high to low). We can do this by introducing the `-r` option.

```Bash
mbge411@login02:~/week01_tutorial/filters$ sort -rnk3 isophya.csv
CKLYU	Pale	2300	0.00	75.81
CKLYU	Pale	2300	0.00	75.81
CKLYU	Pale	2300	0.00	75.81
KALEK	Pale	2100	0.05	76.74
KALEK	Pale	2100	0.05	76.74
KALEK	Pale	2100	0.05	76.74
VRCNK	Pale	2000	0.81	76.92
PIKNK	Pale	1300	0.58	78
IST13	Dark	1000	-0.86	75.66
IST13	Dark	1000	-0.86	75.66
PLVYL	Dark	850	-0.77	73.7
PLVYL	Dark	850	-0.77	73.7
IST06	Dark	450	-0.94	65.9
IST06	Dark	450	-0.94	65.9
```

### wc

`wc`, stands for word count and as the name suggests this command counts words. The default output of this command returns the `number of lines, words and characters` (in that order). If we want to restrict the command to only count `words`, `characters` or `lines` then we can specify this as options.

```Bash
mbge411@login02:~/week01_tutorial/filters$ wc isophya.csv
14 70 373 isophya.csv 
```

Sometimes you want to count only lines `-l`, or words `-w` or characters `-m` or a combination of the two.

```Bash
mbge411@login02:~/week01_tutorial/filters$ wc -l isophya.csv 
14 isophya.csv
```

```Bash
mbge411@login02:~/week01_tutorial/filters$ wc -wm isophya.csv 
70 373 isophya.csv
```

### cut

Usually we work with files that are separated into fields or columns and sometimes we need to filter out a certain field or column to work on just that variable. If this is your need then `cut` is the command for you.

In our sample file we have 5 columns (or fields): `population`, `color`, `altitude`, `body size index` and `humidity`. Let us say we only want to view the `color column`. We can use the following command to achive this.

```Bash
mbge411@login02:~/week01_tutorial/filters$ cut -f2 isophya.csv
Dark
Dark
Dark
Dark
Dark
Dark
Pale
Pale
Pale
Pale
Pale
Pale
Pale
Pale
```

`cut` by default uses the `TAB` character as `field separator`. However, our files might use other separators such as `,` or `space`. For these cases we can specify the separator using the `-d` option. For example, for a file where columns are separated by `spaces` if we wanted to extract `columns 3 and 5` we would use the following command:

```Bash
mbge411@login02:~/week01_tutorial/filters$ cut -d' ' -f3,5 isophya2.csv
450 65.9
450 65.9
850 73.7
850 73.7
1000 75.66
1000 75.66
1300 78
2000 76.92
2100 76.74
2100 76.74
2100 76.74
2300 75.81
2300 75.81
2300 75.81
```

### uniq

`uniq` stands for `unique` and its job is to `remove duplicate lines` from the data. However, it only works if `duplicate lines are adjacent` (i.e. one after the other). Sometimes this is not the case, but this can be fixed by `sorting first` and using `uniq` (we will later see how we can use several commands sequentially via piping).

Let’s say we want to learn how many unique populations we have in our file. For this we simply use the below command.

```Bash
mbge411@login02:~/week01_tutorial/filters$ uniq isophya.csv
IST06	Dark	450	-0.94	65.9
PLVYL	Dark	850	-0.77	73.7
IST13	Dark	1000	-0.86	75.66
PIKNK	Pale	1300	0.58	78
VRCNK	Pale	2000	0.81	76.92
KALEK	Pale	2100	0.05	76.74
CKLYU	Pale	2300	0.00	75.81
```

Like many of the other commands `uniq` also has options. For example, one of the most useful ones is `-c` (count). When we use `uniq` together with the `-c` option, the output will display `how many times each unique item appears` in the file. Let's use the same example as above and see how our output differs.

```Bash
mbge411@login02:~/week01_tutorial/filters$ uniq -c isophya.csv
      2 IST06	Dark	450	-0.94	65.9
      2 PLVYL	Dark	850	-0.77	73.7
      2 IST13	Dark	1000	-0.86	75.66
      1 PIKNK	Pale	1300	0.58	78
      1 VRCNK	Pale	2000	0.81	76.92
      3 KALEK	Pale	2100	0.05	76.74
      3 CKLYU	Pale	2300	0.00	75.81
```

## sed

`sed` is a command that allows us to do a `search and replace` on our data and stands for `Stream Editor`. `sed` is a very powerful and versatile command but here we will only look at its basic function. A basic expression for `search and replace` is like below:

`s/search/replace/g`

The initial `s` stands for `substitute` and designates the action to be performed. Between the `first and second slashes` we place the `pattern` we are `searching`. Between the `second and third slash` we place the `pattern` we wish to `replace` with. `g` stands for `global` and designates we want to replace all instances of the pattern on each line. If we omit the `g` then the command will `replace` the `first instance` of the pattern on each line.

To give an example let us `replace` all instances of `dark` with `blue`:

```Bash
mbge411@login02:~/week01_tutorial/filters$ sed 's/Dark/Blue/g' isophya.csv
IST06	Blue	450	-0.94	65.9
IST06	Blue	450	-0.94	65.9
PLVYL	Blue	850	-0.77	73.7
PLVYL	Blue	850	-0.77	73.7
IST13	Blue	1000	-0.86	75.66
IST13	Blue	1000	-0.86	75.66
PIKNK	Pale	1300	0.58	78
VRCNK	Pale	2000	0.81	76.92
KALEK	Pale	2100	0.05	76.74
KALEK	Pale	2100	0.05	76.74
KALEK	Pale	2100	0.05	76.74
CKLYU	Pale	2300	0.00	75.81
CKLYU	Pale	2300	0.00	75.81
CKLYU	Pale	2300	0.00	75.81
```


## Grep & Regular Expressions

`grep` is a command that will `search` a given file and `print every line` which `contains a given pattern`. In many ways `searching for patterns is the core concept` of bioinformatics so `grep` is one of the, if not the most used, command in a bioinformaticians arsenal. Below are some examples of how grep works.

```Bash
mbge411@login02:~/week01_tutorial/filters$ grep PIKNK isophya.csv
PIKNK	Pale	1300	0.58	78
mbge411@login02:~/ls week01_tutorial/filters$ grep CKLYU isophya.csv 
CKLYU	Pale	2300	0.00	75.81
CKLYU	Pale	2300	0.00	75.81
CKLYU	Pale	2300	0.00	75.81
```

On its own `grep` might not look all that impressive. But when this simple command is combined with the concept of `regular expressions` (re’s) it becomes something very very powerful. So what are `regular expressions`? Regular expressions are similar to `wildcards` and allow us to create complex patterns. A detailed `regular expression` tutorial is out of the scope of this introductory text but there are some great online sources like [this](https://ryanstutorials.net/regular-expressions-tutorial/).

Below is a list of the basic building blocks of `regular expressions`, followed by a few examples of their usage with `grep` so you have some basic idea of how they can be used.

- `. (dot)` -> a single character.
- `?` -> the preceding character matches 0 or 1 times only.
- `*` -> the preceding character matches 0 or more times.
- `+` -> the preceding character matches 1 or more times.
- `{n}` -> the preceding character matches exactly n times.
- `{n,m}` -> the preceding character matches at least n times and not more than m times.
- `[agd]` -> the character is one of those included within the square brackets.
- `[^agd]` -> the character is not one of those included within the square brackets.
- `[c-f]` -> the dash within the square brackets operates as a range. In this case it means either the letters c, d, e or f.
- `( )` -> allows us to group several characters to behave as one.
- `| (pipe symbol)` -> the logical OR operation.
- `g` -> matches the beginning of the line.
- `$` -> matches the end of the line.

Let’s say we want to identify each line that `starts with K`. In the command below we use `^` to designate that the pattern we are searching should be at the `beginning of a line`.

```Bash
mbge411@login02:~/week01_tutorial/filters$ grep '^K' isophya.csv
KALEK	Pale	2100	0.05	76.74
KALEK	Pale	2100	0.05	76.74
KALEK	Pale	2100	0.05	76.74
```

Now let’s do the reverse and search for `lines that end` with the number `2`. Here we will use `$` to designate that the pattern we are searching for should be at the `end of the line`.

```Bash
mbge411@login02:~/week01_tutorial/filters$ grep '2$' isophya.csv
VRCNK	Pale	2000	0.81	76.92
```

How about searching for `any lines` that have the number `4` but the number should be `followed by one or more characters`. Here we use `.` to indicate that there should be `a pattern (any pattern)` directly preceding the number 4 and `+` as a `multiplier`. Note that we also use an argument `-E` to tell `grep` that we will be using `extended regular expressions`. At first it might be confusing to know when or when not to use the `-E` flag when working with `regular expressions` so a useful behavior at first is to always use `-E`.

```Bash
mbge411@login02:~/week01_tutorial/filters$ grep -E '4.+' isophya.csv
IST06	Dark	450	-0.94	65.9
IST06	Dark	450	-0.94	65.9
```

Now let us search for `all occurrences` of the number `6` where it is `repeated twice`. Here we use the multiplier `{2}` to indicate that the character we are searching for should be `repeated twice`.

```Bash
mbge411@login02:~/week01_tutorial/filters$ grep -E '6{2}' isophya.csv
IST13	Dark	1000	-0.86	75.66
IST13	Dark	1000	-0.86	75.66
```

Now let us search for `any line` that has the `pattern IST or KAL`. Here we use `|` to designate `or`.

```Bash
mbge411@login02:~/week01_tutorial/filters$ grep -E 'IST|KAL' isophya.csv
IST06	Dark	450	-0.94	65.9
IST06	Dark	450	-0.94	65.9
IST13	Dark	1000	-0.86	75.66
IST13	Dark	1000	-0.86	75.66
KALEK	Pale	2100	0.05	76.74
KALEK	Pale	2100	0.05	76.74
KALEK	Pale	2100	0.05	76.74
```

Lastly let us search for `all lines` that have a population `beginning` with `letters from J to Z`. 

```Bash
mbge411@login02:~/week01_tutorial/filters$ grep -E '^[J-Z]' isophya.csv 
PLVYL	Dark	850	-0.77	73.7
PLVYL	Dark	850	-0.77	73.7
PIKNK	Pale	1300	0.58	78
VRCNK	Pale	2000	0.81	76.92
KALEK	Pale	2100	0.05	76.74
KALEK	Pale	2100	0.05	76.74
KALEK	Pale	2100	0.05	76.74
```

## Piping & redirecting

Uptill now we have seen a lot of different commands we can use in `Unix/Linux` to manipulate data. However, one of the most powerful features of `Unix/Linux` is that we can join different commands together to do even more powerful stuff. What we essentially do is take the `output` of one command and use it as `input` in another command and we can extend this process until achieving our goal, no matter how many different commands we might need. This process is known as `piping` and coming up with analysis pathways that use multiple commands or programs to manipulate or analyze data are called `pipelines`.

To understand how to build a `pipeline` we must understand a bit about how `Unix/Linux` commands behave. Every command or program we run on the command line has three data streams connected with it.

- `STDIN (0)` -> Standard input (data fed into the program)
- `STDOUT (1)` -> Standard output (data printed by the program, defaults to the terminal)
- `STDERR (2)` -> Standard error (for error messages, also defaults to the terminal)

Once we execute a command or run a program, the output will appear on our screen. This might be suitable for some instances but most of the time we would like to `save` or `write` our output to a separate file. To do this we can use the greater than symbol `>` which takes the `STDOUT` of the command or program and `saves` or `writes` it into a file which we designated.

For example let us `list` the contents of our directory:

```Bash
mbge411@login02:~/week01_tutorial/filters$ ls *
isophya2.csv  isophya.csv  myfile.txt  mylongfile.txt  uarctos.metadata
```

Now instead of listing the contents on the screen let us `write` it to a file called `listoutput`

```Bash
mbge411@login02:~/week01_tutorial/filters$ ls * > listoutput
mbge411@login02:~/week01_tutorial/filters$ cat listoutput 
isophya2.csv
isophya.csv
myfile.txt
mylongfile.txt
uarctos.metadata
```

If we `write` to a file which does not exist, it will be created automatically for us. If we `save` into a file which already exists, however, then it's contents will be `overwritten` by the new output. Here is an example:

```Bash
mbge411@login02:~/week01_tutorial/filters$ wc -l isophya.csv > listoutput 
mbge411@login02:~/week01_tutorial/filters$ cat listoutput 
14 isophya.csv
```

As you can see all of the previous content of `listoutput` has now been `replaced` by the output of the new command (i.e. the number of lines in isophya.csv)

To avoid erasing valuable output we can either write to a separate file each time or we can get the new data to be `appended` (added) onto the same file using the double greater than symbol `>>`.

```Bash
mbge411@login02:~/week01_tutorial/filters$ cat listoutput 
14 isophya.csv
mbge411@login02:~/week01_tutorial/filters$ ls * >> listoutput 
mbge411@login02:~/week01_tutorial/filters$ cat listoutput 
14 isophya.csv
isophya2.csv
isophya.csv
listoutput
myfile.txt
mylongfile.txt
uarctos.metadata
```

So far we have looked at ways to `write` data into files. Now we will look at `piping`, which is a way of `sending data from one command to the next`. For this we will use the symbol `|` which feeds the `output` of the command on the `left` as `input` to the command on the `right`.

In the example below we will `list` all the contents of our directory using the command `ls` and then display only the first 3 contents using the command `head`.

```Bash
mbge411@login02:~/week01_tutorial/filters$ ls * | head -3 
isophya2.csv
isophya.csv
listoutput
```

We can extend this as long as we want. For example, we can list all contents our directory, then get the first 3 and then get the last item on the reduced list of three items:

```Bash
mbge411@login02:~/week01_tutorial/filters$ ls * | head -3 | tail -1 
listoutput
```

`Piping` is an effective tool for getting quick answers about our data. For example, let's say I want to know how many `Pale` individuals are in the file `isophya.csv`. I can answer this question with a simple command like below.

```Bash
mbge411@login02:~/week01_tutorial/filters$ grep Pale isophya.csv | wc -l 
8
```

If I want to `write` the answer to a `file` instead of on the screen, then I can add an additional step and save it to a file:

```Bash
mbge411@login02:~/week01_tutorial/filters$ grep Pale isophya.csv | wc -l > count.pale
mbge411@login02:~/week01_tutorial/filters$ cat count.pale 
8
```

## awk

`awk` is a `comprehensive programming language` for text-processing on `Unix/Linux` and is the cornerstone of `Unix/Linux` shell programming. Unfortunately we do not have time to go into all of it here. The things you can do with `awk` are limitless and I will only give a few example to get you stared and leave it up to you to discover the rest. An introductory tutorial can be found [here](https://likegeeks.com/awk-command/).

On of the most basic thing we can do with awk is to print a specific column or columns. For example we may wish to print only the 1st and 3rd columns in our isophya.csv file.

```Bash
mbge411@login03:~/course_content/week01_tutorial/filters$ awk '{print $1, $3}' isophya.csv
IST06 450
IST06 450
PLVYL 850
PLVYL 850
IST13 1000
IST13 1000
PIKNK 1300
VRCNK 2000
KALEK 2100
KALEK 2100
KALEK 2100
CKLYU 2300
CKLYU 2300
CKLYU 2300
```

We can also print in any order we want

```Bash
mbge411@login03:~/course_content/week01_tutorial/filters$ awk '{print $3, $1}' isophya.csv
450 IST06
450 IST06
850 PLVYL
850 PLVYL
1000 IST13
1000 IST13
1300 PIKNK
2000 VRCNK
2100 KALEK
2100 KALEK
2100 KALEK
2300 CKLYU
2300 CKLYU
2300 CKLYU
```

We can also do simple calculations between columns using awk. For example let us say we want to divide the 3rd column in isophya.csv with the 5th column and print this as a new variable. 

```Bash
mbge411@login03:~/course_content/week01_tutorial/filters$ awk 'c=$5/$3 {print c}' isophya.csv
0.146444
0.146444
0.0867059
0.0867059
0.07566
0.07566
0.06
0.03846
0.0365429
0.0365429
0.0365429
0.0329609
0.0329609
0.0329609
```

awk can also be used to filter our dataframe according to certain values. For example we might want to only work with populations that are found above 1000 meters and in environments where humidity is lower than 76%.

```Bash
mbge411@login03:~/course_content/week01_tutorial/filters$ awk '$3>1000 && $5<76' isophya.csv
CKLYU	Pale	2300	0.00	75.81
CKLYU	Pale	2300	0.00	75.81
CKLYU	Pale	2300	0.00	75.81
```


## THE END

Congratulations you have reached the end of your first `Unix/Linux` tutorial!! Be warned that you have only scratched the surface of what you can do with this very powerful system. I have intentionally left some things out  because a good portion of learning the `Unix/Linux` and the `Bash` environment is going `online` and finding solutions to your problems. There are many online sources like [Stackoverflow](https://stackoverflow.com/), the Unix and Linux [Forums](https://www.unix.com/), or [explainshell](https://explainshell.com/) and you should make a habit of checking them out regularly. 

A good friend once told me `to be a good coder you only need to know how to ask google the right questions`.

Now that you have reached the end of this tutorial, see if you can complete some of the challenges in this week's homework. Be warned not all answers have been covered here, so happy hunting!!!







