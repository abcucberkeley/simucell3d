# Use a base image with necessary dependencies
FROM ubuntu:latest
FROM gcc:11.4
FROM python:3.10
ENV RUNNING_IN_DOCKER=TRUE

# Make bash colorful https://www.baeldung.com/linux/docker-container-colored-bash-output   https://ss64.com/nt/syntax-ansi.html
ENV TERM=xterm-256color
RUN echo "PS1='\e[97m\u\e[0m@\e[94m\h\e[0m:\e[35m\w\e[0m# '" >> /root/.bashrc

# Set working directory
WORKDIR /app

# Install CMake
RUN apt-get update && apt-get install -y cmake

# Install requirements. Don't "apt-get upgrade" or else all the NVIDIA tools and drivers will update.
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
  sudo \
  htop \
  cifs-utils \
  winbind \
  smbclient \
  sshfs \
  iputils-ping \
  && rm -rf /var/lib/apt/lists/*

# Git-lfs install
RUN curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash && apt-get install git-lfs && rm -rf /var/lib/apt/lists/*



# Create a build directory and run CMake to configure the project
RUN mkdir -p container_build 

#Create a foder to store the simulation results
RUN mkdir -p simulation_results

# Copy source code into the container (and avoid copying in all the simulation results)
COPY src src
COPY lib lib
COPY test test
COPY include include
COPY CMakeLists.txt .
COPY main.cpp .
COPY scripts scripts

# Set working directory to the build directory
WORKDIR /app/container_build
# Prepare the build
RUN cmake -DENABLE_PYTHON_BINDINGS=TRUE -DPYTHON_EXECUTABLE=$(which python3) -DCMAKE_BUILD_TYPE=Release ..

# Build the project
RUN make -j4


# Add python packages
RUN pip install --upgrade pip
RUN pip install --no-cache-dir --progress-bar off pyvista tqdm ipython line_profiler_pycharm pytest -r /app/scripts/python_bindings/requirements.txt &&  pip cache purge || true

ARG USERNAME=vscode
ENV USER=${USERNAME}
ARG USER_UID=1000
ARG USER_GID=1000

# Create the user
RUN   groupadd --gid $USER_GID $USERNAME && \
    groupadd --gid 1001 vscode_secondary && \
    useradd -l --uid $USER_UID --gid $USER_GID -G 1001 -m $USERNAME && \
    #
    # [Optional] Add sudo support. Omit if you don't need to install software after connecting.
    echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME || true

# [Optional] Set the default user. Omit if you want to keep the default as root.
# USER $USERNAME


# Define the command to run when the container starts
ENTRYPOINT ["/app/container_build/simucell3d"]