# Some git tricks

### Checkout a remote branch

```
git fetch
git branch -v -a    # to visualise the branches
git switch -c <remote_branch_name> origin/<remote_branch_name>
```

### Merge a local branch with another branch

Situation: I have a local branch (e.g. `anelastic`) that I want to modify without damaging other branches. There exists a `origin/anelastic` to which I push regularly. But I need updates being developed in another branch (e.g. `origin/feature/ito`).

```
git fetch; git merge origin/feature/ito
```
