# Gamma Ray Detectors in Space

For this experiment, you will have to build a pipeline (or, possibly, a small
number of related but distinct pipelines) with which to process the spectral
data that you have collected.

We will use GitHub to help monitor your progress. You will maintain your own
code in this repository, keeping the repository up-to-date with your latest
work. In the lab sessions, we will discuss your code with you, offering hints
and suggestions for improvements and development.

## Workflow

You should create a local copy of your repository on your machine. You can do
this from the command prompt by using the command

```
git clone https://github.com/UCDSatelliteSubsystems/...
```
where you will replace the ellipsis with the rest of the URL path to your
repository. Alternatively, the GitHub desktop app will help you to do this.

This local clone of the repository is completely independent of the remote
repo that you have set up on GitHub. If you make changes, they will only
affect your local repo until you share them with the remote server.


### Push & Pull

The `git push` and `git pull` commands allow you to synchronise your local
repo with the remote repo.

* You use `git push` to push changes you have made locally on to the remote
  server.
* You use `git pull` to pull your local repo up to the latest
  version on the remote server.

This is usually a simple task, but occasionally you will find that your local
copy and the remote copy are in conflict, usually because commits have been
made to both that make different changes to the same piece of code. If this
happens, Git will ask you to resolve these conflicts, and then commit and push
the newly synchronised version.

The best way to avoid conflicts is to `push` and `pull` often. You should pull
the latest version at the start of each development session, and commit and
push your changes just before you finish up.

### Branches

Git uses branches to allow you to maintain different versions of a codebase
simultaneously. For example, a software devloper might have one branch
containing the latest stable release version, and another branch containing
the bleeding-edge development version, which may have new features that are
under-construction.

You should maintain at least two branches in your repository:

* `main` should contain the latest version that you have submitted. This
  should be a fairly stable branch, perhaps changing only once a week.
* `dev` should be the development branch, where you do your work and make
  commits as you make progress.

You can list the branches in your repo by typing `git branch -a`.

You can create this new development branch in your local repo with the command
`git checkout -b dev`. Of course, you can create other branches by giving them
different names.

You will then need to push this branch to the remote server, so we can keep
track of it. You will need to tell the remote server that there is a new
branch, which you can do using:

```
git push --set-upstream origin dev
```

You can merge branches together locally. You might need to do this if you
create further branches off the `dev` branch. However, within this workflow,
you should not locally merge `dev` and `main`. Instead, you will use pull
requests.


### Pull Requests

Before each of the review dates below, you should open a pull request to merge
your development branch back into the main branch. You can do this using the
GitHub web interface. On the repo page, select Pull Request and then New Pull
Request.

When you create the pull request, you will be asked to explain the content of
the request. Briefly summarise the progress made and the next steps. We will
use this to guide our discussion in the lab.

## Timeline

As an approximate timeline,


| Date    | Acheivement   |
| ------- | --------------- |
| 1 / 10  | Read in and plot spectra | 
| 8 / 10  | Fit models to peaks | 
| 15 / 10 | Calculate parameters (eg., resolution, efficiency) | 
| 22 / 10 | Completed analysis pipeline  |



## Some Dos and Don'ts

* Do give useful commit messages.

* Don't upload your raw data or your figures. The repository is for the code
  you are writing. It is not for the data you are applying it to. If you mix
  the two, then (a) the repository will get quite large, quite quickly, and
  (b) the repository will get messy and it will be harder to keep track of the
  important things.