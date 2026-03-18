import json
import logging
import os
import re
import subprocess
import sys
from pathlib import Path

# Constants for default values
DEFAULT_PROD_IMAGE_BASE = 'australia-southeast1-docker.pkg.dev/cpg-common/images'
DEFAULT_DEV_IMAGE_BASE = 'australia-southeast1-docker.pkg.dev/cpg-common/images-dev'
DOCKERFILE_NAME = 'Dockerfile'
CONTAINER_NAME = 'popgen_genotyping'

# Regex for extracting version from Dockerfile: ENV VERSION=1.0.0
VERSION_PATTERN = re.compile(r'^\s*ENV\s+VERSION\s*=\s*([^\s]+)', re.MULTILINE)

logging.basicConfig(level=logging.ERROR)


def extract_version_from_file(file_path: str | Path) -> str | None:
    """
    Extract the version from a Dockerfile by searching for a line like:
    ENV VERSION=1.0.0
    """
    path = Path(file_path)
    if not path.exists():
        return None

    content = path.read_text()
    match = VERSION_PATTERN.search(content)
    return match.group(1) if match else None


def list_image_tags(full_image_name: str) -> list[dict]:
    """
    Query GCP to list tags for the given image name.
    """
    cmd = [
        'gcloud',
        'container',
        'images',
        'list-tags',
        full_image_name,
        '--format=json',
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)  # noqa: S603
        return json.loads(result.stdout)
    except (subprocess.CalledProcessError, json.JSONDecodeError) as e:
        logging.error(f'Failed to list tags for {full_image_name}: {e}')
        return []


def get_next_version_tag(container_name: str, version: str) -> str:
    """
    Determines the next available version suffix for the extracted version
    by querying both production and development image repositories.
    """
    base_prod = os.environ.get('GCP_BASE_IMAGE', DEFAULT_PROD_IMAGE_BASE)
    base_dev = os.environ.get('GCP_BASE_ARCHIVE_IMAGE', DEFAULT_DEV_IMAGE_BASE)

    tags_list = []
    for base_path in [base_prod, base_dev]:
        full_name = f'{base_path}/{container_name}'
        tags_list.extend(list_image_tags(full_name))

    # Pattern to match tags like '1.0.0-5'
    tag_pattern = re.compile(rf'^{re.escape(version)}-(\d+)$')

    max_suffix = 0
    for entry in tags_list:
        for tag in entry.get('tags', []):
            match = tag_pattern.match(tag)
            if match:
                max_suffix = max(max_suffix, int(match.group(1)))

    return f'{version}-{max_suffix + 1}'


def main():
    current_version = extract_version_from_file(DOCKERFILE_NAME)
    if current_version is None:
        raise NotImplementedError('The Dockerfile needs to contain a version string in the format "ENV VERSION=x.x.x"')

    # Determine the next available tag based on current_version.
    new_tag = get_next_version_tag(CONTAINER_NAME, current_version)

    matrix = {'include': [{'name': CONTAINER_NAME, 'tag': new_tag}]}

    # Output for CI consumption
    json_output = json.dumps(matrix, separators=(',', ':'))
    print(json_output, file=sys.stderr)
    print(json_output, end='')


if __name__ == '__main__':
    main()
