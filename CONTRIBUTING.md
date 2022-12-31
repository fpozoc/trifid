<h1>Contributor Guidelines</h1>

- [Bug Reports](#bug-reports)
- [Feature Requests](#feature-requests)
- [Simple Issue template](#simple-issue-template)
- [Documentation Contributions](#documentation-contributions)
- [Code Contributions](#code-contributions)
- [Submitting changes](#submitting-changes)
  - [Git branches](#git-branches)
  - [Git commit messages](#git-commit-messages)

This document lays out guidelines and advice for contributing to this project.
If you're thinking of contributing, please start by reading this document and getting a feel for how contributing to this project works.

The guide is split into sections based on the type of contribution you're thinking of making.

## Bug Reports

Bug reports are hugely important! Before you raise one, though, please check through the [`GitHub issues`](https://github.com/fpozoc/trifid/issues), **both open and closed**, to confirm that the bug hasn't been reported before.

When filing an issue, make sure to answer these questions:

- Which Python version are you using?
- Which version of `trifid` are you using?
- What did you expect to see?
- What did you actually see?
- What are you going to do?

The best way to get your bug fixed is to provide a test case, and/or steps to reproduce the issue.

## Feature Requests

If you believe there is a feature missing, feel free to raise a [`Feature request`](https://github.com/fpozoc/trifid/issues/new?assignees=&labels=&template=feature_request.md&title=) to answer:

- A clear and concise description of what the problem is.
- The solution you would like.
- Alternatives you have considered.
- Additional context.

## Simple Issue template

If you don't feel confortable with the previous options, you can also use the [`Simple Issue template`](https://github.com/fpozoc/trifid/issues/new?assignees=&labels=&template=simple-issue-template.md&title=)

## Documentation Contributions

Documentation improvements are always welcome!

You do not have to setup a development environment to make small changes to the docs. Instead, you can `edit files directly on GitHub`_ and suggest changes.

## Code Contributions

If you intend to contribute code, do not feel the need to sit on your contribution until it is perfectly polished and complete. It helps everyone involved for you to seek feedback as early as you possibly can. Submitting an early, unfinished version of your contribution for feedback can save you from putting a lot of work into a contribution that is not suitable for the project.

## Submitting changes

1. Create an issue.
2. Assign the issue to a developer and create a branch.
3. Create a Pull Request.
4. The admin (first creator) of the repository will review the Pull Request.

### Git branches

- Always perform work in a feature branch.
- It is better to branch out from `develop` branch.
- Delete local and remote feature branches after merging.

### Git commit messages

- Separate the subject from the body with a newline between the two (if the body exists).
- Limit the subject line to 50 characters and Wrap the body at 72 characters.
- Capitalize the subject line.
- Do not end the subject line with a period.
- Use imperative mood in the subject line.

Example:

```sh
git commit -m "Update getting started documentation"
```

(If applied, this commit will **update getting started documentation**)
