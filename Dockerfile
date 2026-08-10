FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_flow:1.3.1

# CI's get_version.py reads this line to derive the image tag (VERSION-<n>);
# it must stay a literal semver, so the VERSION build-arg CI passes is unused.
ENV VERSION=0.1.0

# Set the working directory
WORKDIR /app

# Copy the build configuration and README
COPY pyproject.toml README.md LICENSE ./

# Copy the source code
COPY src/ ./src/

# Install the project and its dependencies
# Using --no-cache-dir to keep the image small
RUN pip install --no-cache-dir .

# The entrypoint can be set to the run_workflow script if desired
# ENTRYPOINT ["run_workflow"]
