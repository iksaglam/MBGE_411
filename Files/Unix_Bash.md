
# Introduction to Linux and Bash 
___

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
