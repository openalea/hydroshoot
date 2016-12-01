============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given.


You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at issues_.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug"
is open to whoever wants to implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with "feature"
is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

pkglts could always use more documentation, whether as part of the
official pkglts docs, in docstrings, or even on the web in blog posts,
articles, and such.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at issues_.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Get Started!
------------

Ready to contribute? Here's how to set up `hydroshoot` for local development.

1. Fork the `hydroshoot` repo on GitHub.
2. Clone your fork locally::

    $ git clone git@github.com:your_name_here/hydroshoot.git

3. Install your local copy into a virtualenv. Assuming you have virtualenv_
installed, this is how you set up your fork for local development::

    $ virtualenv dvlpt
    $ dvlpt/script/activate
    (dvlpt)$ python setup.py develop

4. Create a branch for local development (wip stands for work in progress)::

    (dvlpt)$ git checkout -b wip_name-of-your-bugfix-or-feature

   Now you can make your changes locally.

5. When you're done making changes, check that your changes pass flake8 and the
tests, including testing other Python versions with tox::

    (dvlpt)$ cd hydroshoot
    (dvlpt) hydroshoot$ flake8
    (dvlpt) hydroshoot$ nosetests
    (dvlpt) hydroshoot$ tox

   To get flake8 and tox, just pip install them into your virtualenv.

6. Commit your changes and push your branch to GitHub::

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin wip_name-of-your-bugfix-or-feature

7. Submit a pull request through the GitHub website.

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add the
   feature to the list in README.rst.
3. The pull request should work for Python 2.7, 3.4, and 3.5. Check
   https://travis-ci.org/moi/hydroshoot/pull_requests
   and make sure that the tests pass for all supported Python versions.

Tips
----

To run a subset of tests::

    $ nosetests test/test_XXX



.. _issues: https://github.com/moi/hydroshoot/issues
.. _virtualenv: https://pypi.python.org/pypi/virtualenv
