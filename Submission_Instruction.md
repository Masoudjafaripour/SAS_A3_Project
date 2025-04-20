# Submission Instruction

## Login with your GitHub account

Login to the [competition website](https://gppc.search-conference.org/grid) with a GitHub account, and we will automatically create a private GitHub submission repo for you.
The repo will be the place that you submit codes to. You can click "My Repo" to open your GitHub submission repo page.

## Clone your submission repo

Clone your submission repo to your local machine. The repo contains starter code to help you prepare your submission.

```
$ git clone git@github.com:your_submission_repo_address
$ cd your_submission_repo
```

## Implement your algorithm

Read the Problem_Definition.md to check the problem definitions.

Read the start kit README.md to know what files you should (not) modify, and where you should implement your algorithm.

When your implementation is ready to evaluate, add, commit, and push all your changes by running the following commands in your local repo:
```
$ git add *
$ git commit -m "Some message to describe this commit."
$ git push origin
```

## Evaluate your algorithm

Once you commit your implementation to your submission repo you can begin the evaluation process. From the competition website, navigate to the any-angle track. Here you will be able to see details of your entry and submission history (make sure you are logged in!). 

If your implementation support multi-thread preprocessing, tick the "My Entry Support Multi-thread Preprocessing" option.
If ticked, the preprocessing server will assign 4 CPUs for preprocessing (otherwise only 1 CPU is assigned).

You can evaulate a specific branch under "Evaluate the Branch", this allows for different algorithms or variant features.

Then click the "Evaluate my codes" button on the competition website to evaluate your new submission.

## Track evaluation progress and history submission

Click the "My Submissions" button to see your submission history.
Click an entry in the history to see details and track the progress of a running evaluation.
