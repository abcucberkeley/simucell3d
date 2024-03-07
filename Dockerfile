# Use a base image with necessary dependencies
FROM ubuntu:latest
FROM gcc:11.4
FROM python:3.10


# Set working directory
WORKDIR /app

# Copy source code into the container
COPY . .

# Install CMake
RUN apt-get update && apt-get install -y cmake

# Create a build directory and run CMake to configure the project
RUN mkdir -p container_build 

#Create a foder to store the simulation results
RUN mkdir -p simulation_results

# Set working directory to the build directory
WORKDIR /app/container_build

# Prepare the build
RUN  cmake -DENABLE_PYTHON_BINDINGS=TRUE -DPYTHON_EXECUTABLE=$(which python3) -DCMAKE_BUILD_TYPE=Release .. 

# Build the project
RUN make -j4

# Define the command to run when the container starts
ENTRYPOINT ["/app/container_build/simucell3d"]