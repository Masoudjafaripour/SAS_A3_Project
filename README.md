# Grid-Based Path Finding Competition Starter Kit

## Grid Track

Website [https://gppc.search-conference.org/grid](https://gppc.search-conference.org/grid)

## Submission Instruction

Read the [Submission_Instruction.md](https://github.com/gppc-dev/startkit-anyangle/blob/master/Submission_Instruction.md)

## Problem Definition

Read the [Problem_Definition.md](https://github.com/gppc-dev/startkit-classic/blob/master/Problem_Definition.md)

# TLDR

* Participants push commits to their repositories and the server will pull, compile, run and evaluate the current head of the main branch.

  * participants can add new repo via our web interface
  * participants must specify their dependency in `apt.txt` (we provide a sample in `startkit`)
  * server will build a docker image to compile, run and evaluate submissions

# Start Kit

For those who using c++ (most participants), your submission must include the following files:

| File name             | Description                                                     | Modifiable |
| --------------------- | --------------------------------------------------------------- | ---------- |
| `main.cpp`            | Define API of executable, compiled to `run`                     | no         |
| `Timer.h`             | Define timer                                                    | no         |
| `Timer.cpp`           | Define timer                                                    | no         |
| `ScenarioLoader.h`    | GPPC scenario file parser & loader                              | no         |
| `ScenarioLoader.cpp`  | GPPC scenario file parser & loader                              | no         |
| `GPPC.h`              | Commonly used code for GPPC                                     | no         |
| `Entry.h`             | Define functions prototypes that will be used by `main.cpp`     | no         |
| `validator/*`         | Validator code                                                  | no         |
| `version.txt`         | Startkit version                                                | no         |
| `Entry.cpp`           | Implementations                                                 | yes        |
| `compile.sh`          | Compile your code to executable `run`                           | yes        |
| `apt.txt`             | Define dependency, will be used by server to build docker image | yes        |

Your program may generated the following files:

| File name     | Description                                                                  | Optional |
| ------------- | ---------------------------------------------------------------------------- | -------- |
| `run`         | compiled executable, will be called in evaluation                            | no       |
| `run.stdout`  | stdout is redirected to here. Computed paths goes to stdout for validation.  | no       |
| `run.stderr`  | stderr redirected to here                                                    | no       |
| `run.info`    | stores some run time information                                             | no       |
| `result.csv`  | stores query information, including time cost, path length, etc.             | no       |
| `index_data/` | if your algorithm has pre-computation, all produced data must be here        | yes      |

We provide `A*` in c++ as an example.

For those who using other languages, you may not need to include `*.cpp` and `*.h` files in the above table, but others are still required.
Notice that `run` is an untracked file by default (see in `.gitignore`), if you put code in this file, make sure you also modify the `.gitignore`.

## Your Implementation
* Implement `PreprocessMap`, `PrepareForSearch`, and `GetPath` functions in `Entry.cpp`. See examples and detailed documentations in `Entry.cpp`.
* Specify your dependency packages in `apt.txt`. The packages must be available for installation through `apt-get` on Ubuntu 22.
* Modify `compile.sh` and make sure your code can be compiled by executing this script.

## Run the Program
* `./run -pre <map> none` Run in preprocessing mode. The program should preprocess the given map and store the preprocessing data under `index_data/`.
* `./run -check <map> <scen>` Run in validation mode. The output will be validated. Each entry of the `run.stdout` will be marked as `valid` or `invalid-i`, where `i` indicate which segment of the path is invalid.
* `./run -run <map> <scen>` Run in benchmark mode. The benchmark results are written to `result.csv`.

## Customise Program Runtime

Environmental variables are defined to enable features not strictly required for development.
They are listed below:

* `GPPC_REDIRECT_OUTPUT`: redirects `stdout`/`stderr` to files, as detailed in I/O Setup section.
* `GPPC_MEMORY_TRACK`: prints memory usage into `run.info` file, available on Linux only.

# Details on the server side

For those who **want to build local testing workflow** or **not using c/c++**, this section might be helpful.

## I/O Setup

* All `stdout` from program are redirected to a file `run.stdout`. Your codes should not print any debug info to `stdout`.  

* All `stderr` are redirected to a file `run.stderr`

* The results of benchmark (i.e. `./run -run <map> <scen>`) are written to `result.csv`

* All these files are in docker, and will backup to server so that we can hide/reveal information to participants.


## Execution Setup

* The server run all scripts under certain constraints, e.g. time limit, memory limit.

* A execution will be killed by server when it exceeds limits, and participants will get `Time Limit Exceed` / `Memory Limit Exceed` verdict.

* If an execution crashed, participants will get `Run Time Error` verdict with few or none hints.

* Participants can submit their solution and execute in `dev` mode for testing purpose, in `dev` mode:
  * we only test a small set of simple map
  * we will reveal all `stdout` and `stderr` to participants
  * this mode is to eliminate compile errors and illegal output format.
  * although participants can do this locally, we encourage everyone to submit on `dev` mode at least once.


## Evaluation Workflow

### Overview
1. Build docker image based on the apt.txt in submission repo.
2. Start the container in background.
3. Run pre-processing for debug maps.
4. Run validation for debug scenarios.
5. Run pre-processing for benchmark maps.
6. Run validation for benchmark scenarios.
7. Run benchmark for benchmark scenarios.
8. Submit final result.


# Test Your Implementation in a Docker Environment

* Install latest docker release on your machine, [https://docs.docker.com/engine/install/](https://docs.docker.com/engine/install/).

* Setting up environment using `RunInDocker.sh`:
  * In the root of your code base, run command `./RunInDocker.sh`. This script will automatically generate a Dockerfile to build the docker image.
  * It will copy your codes to the Docker Environment, install dependencies listed in `apt.txt` using apt-get, and compile your code using `compile.sh`.
  * You are inside the docker container when the script finishes.
  * You can run the compiled program in side docker container now.
  * The docker image name `<image name>` is `gppc_image` and the container name `<container name>` is `gppc_test`.
  * Exit the container with command: `exit`.
* Or you also can set up environment manually:
  * Prepare your Dockerfile, which uses `ubuntu:jammy` as the base image.
  * The docker image working directory should be set to the directory where executables are. The docker image should copy user implementation into the image.
  * `docker build -t <image name> <dockerfile>`: build docker image based on a default/user-defined Dockerfile.
  * Building executable:
    * The Dockerfile include a RUN command that running a user provided compile script `./compile.sh` to build executable.
  * `docker run -it --name <container name> <image name>`: start the container interactively.

* Start an existing container:
  * In background: `docker container start <container name>`
  * Interactively: `docker container start -i <container name>`

* The docker container is started in background, you can run commands from the outside of the docker container (treat docker container as executable).
 
  * Use prefix: `docker container exec <container name>` for all following commands
 
  * Running preprocessing:
    * `<prefix> ./run -pre ${map_path} none`: run user provided preprocessing codes.

  * Running executable with the validator
    * `<prefix> ./run -check ${map_path} ${scenario_path} `
    * With `-check` flag, the output will be validated by a validator. Each entry of the output will be marked as `valid` or `invalid-i`, where `i` indicate which segment of the path is invalid.

  * Running executable for benchmarking: `<prefix> ./run -run ${map_path} ${scenario_path}`
    * we will track time/memory usage and publish results
 
  * All outputs are stored inside the container. You could copy files from docker container. For example: `docker cp gppc_test:/GPPC2021/codes/run.stdout ./run.stdout`, which copies `run.stdout` to your current working directory.
