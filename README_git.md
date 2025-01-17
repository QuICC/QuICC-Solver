# Some git tricks

### Checkout a remote branch

```
git fetch
git branch -v -a    # to visualise the branches
git switch -c <remote_branch_name> origin/<remote_branch_name>
```

### Merge branches

Situation: I have a local branch (e.g. `anelastic`) that I want to modify without damaging other branches. There exists a `origin/anelastic` to which I push regularly. But I need updates being developed in another branch (e.g. `origin/feature/ito`).

```
git fetch; git merge origin/feature/ito
```

Merge also works to merge your work on local branches. Imagine I am on `main` (on my local machine) and want to implement a new feature. This would be the workflow:

```
git checkout -b new_feature
# work on the new branch and implement the new feature
git add <path_to_add>
git commit -m "<commit_message>"
git checkout main           # go back on main
git merge new_feature       # merge my feature in main
```
