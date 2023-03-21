# ONSAS.jl contributor guide

Welcome to ONSAS.jl contributor brief guide. Everyone is kindly welcome and of  course there are several things to do. Contributing into this does not require Julia experience or even advanced programming knowledge.

First we strongly recommend to **join the regular community meets or [via slack][slack-link]**. Also you can mail us at: [mvanzulli@fing.edy.uy][mailto-mvanzulli], [mforets@gmail.com][mailto-mforets] and [jorgepz@fing.edu.uy][mailto-jorgepz]

Here you can find things about:

 - [Reporting issues](#issues)
 - [Docs](#docs)
 - [Code changes](#code-changes)
 - [Authors](#authors)

If you want to contribute to a software package, but you're new to open source development you can check on [this page][contributing] to get started. Additionally, if you're interested in contributing to a Julia package, you may find the video [Open source, Julia packages, git, and GitHub][tim-git] to be a valuable resource. 

## Issues

If you have found a bug or a problem with ONSAS.jl you can open an [issue][new-issue]. Try to include as much information about the problem as possible and preferably some code that can be copy-pasted to reproduce it (see [How to create a Minimal, Reproducible Example][create-rep-example]). If you can identify a fix for the bug you can submit a pull request without first opening an issue, see [Code changes](#code-changes).

## Docs

Starting to contribute to a new project's documentation is an excellent way to begin. As a new user, you bring a fresh perspective and can easily identify what areas need better documentation or explanation. If something is confusing to you, chances are other users have similar questions. It's important to remember that even small contributions, such as correcting typos, are valuable. If you're not sure where to begin, you can check out the [open issues][open-issues] for specific areas that need improvement.

Making small changes is simple using GitHub's web interface (refer to [Editing files][gh-edit-files]). For every page in the documentation, there is an Edit on GitHub button at the top, which redirects you to the relevant source file. If you need guidance on how to make Julia documentation better, check out the step by step tutorial: [Making Julia documentation better][tim-doc]. 

Making larger changes is easier locally after cloning the repository. For this process the documenting [Julia manual][julia-doc] and [documentation for `Documenter.jl`][documenter] are essential. To preview the docs you should move into the `./docs/` folder and run: 

```bash
julia --project=.
```

Followed by:

```julia
include("make.jl")
```

A new `docs/build/` folder will be automatically created, to preview the docs open the `index.html` file. 


## Code changes

In case you encounter a bug or issue with ONSAS.jl, you can create an [issue][open-issues] by providing as much information as possible about the problem. Ideally, include some code that can be easily copied and pasted to replicate the issue (refer to [How to create a Minimal, Reproducible Example][create-rep-example]). If you're able to identify a solution for the bug, please create an [issue][open-issues] and then a [pull request][open-pull-request] solving it.  


## Authors 

The authorship of the code is based on the criteria defined by the [JOSS journal][joss]. The co-authors have collaborated in tasks such as: design, important new features development or extensive documentation contributions.


[documenter]: https://juliadocs.github.io/Documenter.jl/
[tim-git]: https://youtu.be/cquJ9kPkwR8
[tim-doc]: https://youtu.be/ZpH1ry8qqfw
[gh-edit-files]: https://docs.github.com/en/repositories/working-with-files/managing-files/editing-files#editing-files-in-another-users-repository
[contributing]: https://contributing.md/
[open-issues]: https://github.com/ONSAS/ONSAS.jl/issues/new
[open-pull-request]: https://github.com/ONSAS/ONSAS.jl/compare
[create-rep-example]: https://stackoverflow.com/help/minimal-reproducible-example
[julia-doc]: https://docs.julialang.org/en/v1/manual/documentation/
[mailto-jorgepz]: mailto:jorgepz@fing.edu.uy
[mailto-mvanzulli]: mailto:mvanzulli@fing.edu.uy
[mailto-mforets]: mailto:mforets@gmail.com
[slack-link]: https://app.slack.com/client/T04QWNG5T2Q/C04R6TMDV0R
[joss]: https://joss.theoj.org/