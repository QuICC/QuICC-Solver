Get the code {#pManGet}
============

The sources are under Git to provide versioning. If you're not familiar with Git, the reference book [Pro Git](http://git-scm.com/book/) is available online.

1. Send me your public SSH key(s)

2. Get the sources:

   #$>git clone gitolite@base13.ch:QuICC.git

3. Go into new diretory

   #$>cd QuICC

4. Get external modules (currently only Eigen):

   #$>git submodule init

   #$>git submodule update

Other useful stuff
------------------

- Commit your part of your changes

   #$>git add_modified_files

   #$>git commit

- Commit all your changes (ignores new files)

   #$>git commit -a

- Update your working directory

   #$>git pull

- Start a new remote branch:

   1. Create new branch:

      #$>git push origin origin:refs/head/new_branch_name

   2. Make sure everything is up-to-date:

      #$>git fetch origin

   3. Check what you have done:

      #$>git branch -r

   4. Track new branch:

      #$>git checkout \-\-track -b bew_branch_name origin/new_branch_name

   For more details see [Start a new remove branch](http://www.zorched.net/2008/04/14/start-a-new-branch-on-your-remote-git-repository/).
