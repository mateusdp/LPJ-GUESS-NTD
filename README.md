# LPJ-GUESS for Senckenberg

## About This Git Repository

This Git repository is for internal use in the [“BIMODAL” working group at Senckenberg BiK-F][bimodal].
It is forked from and regularly synced with the `trunk` of the LPJ-GUESS Subversion repository by the 

Some Git branches are prefixed with `svn-`.
They “mirror” the corresponding branch of the Subversion repository in Lund.
For example, the file tree of the Git branch `svn-trunk` is exactly equivalent to the Subversion `trunk` in its respective revision.
See [this tutorial on “Git+SVN”][git+svn] to learn more about the concept of using Git on top of a Subversion file tree.

The `master` and `develop` branches are protected and can only writable for users with “Maintainer” permission.
See [here][main-branches-explained] to learn about the basics of main branches and feature branches in Git.

Changes to the `master` branch are tracked in the [CHANGELOG.md](CHANGELOG.md), which follows the [keep-a-changelog][] standard.

## About LPJ-GUESS
In the [reference](reference/) folder you will find resources to familiarize yourself with the code base, how to use it, and the scientific background of the model.

You can compile and browse the [Doxygen](https://doxygen.nl) documentation by executing `doxygen doxygen.conf` within the `doxygen/` directory.
This should generate the file `doxygen/output/html/indx.html`, which you can then open with your web browser.

You will find more resources on the [LPJ-GUESS website][] and on the [LPJ-GUESS Wiki][].

## Repository Structure
The LPJ-GUESS directory structure is explained in the original [readme.txt](readme.txt).

## License
The LPJ-GUESS code does not have a license and is thus implicitly proprietary.

The licensing of code contributions from the BIMODAL group is as of yet undefined and therefore also implicitly under exclusive copyright by the contributors or their employers.

[This text][no-license] explains the implications of there being no license for this project.



[bimodal]: https://www.senckenberg.de/en/institutes/sbik-f/quantitative-biogeography/
[git+svn]: https://lostechies.com/derickbailey/2010/02/03/branch-per-feature-how-i-manage-subversion-with-git-branches/
[keep-a-changelog]: https://keepachangelog.com/en/1.0.0/
[LPJ-GUESS website]: http://web.nateko.lu.se/lpj-guess/index.html
[LPJ-GUESS Wiki]: http://stormbringer.nateko.lu.se/public/simba_mirror/Wiki/LPJ_GUESS
[main-branches-explained]: http://web.nateko.lu.se/lpj-guess/index.html
[no-license]: https://choosealicense.com/no-permission/