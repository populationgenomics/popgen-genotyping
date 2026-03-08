"""
Job logic for merging individual BCFs into a cohort PLINK2 dataset.
"""

import json
from pathlib import Path
from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch

if TYPE_CHECKING:
    from hailtop.batch.job import Job


def run_cohort_bcf_to_plink(
    bcf_paths: dict[str, str],
    output_prefix: str,
    job_name: str = 'cohort_bcf_to_plink',
) -> 'Job':
    """
    Run the PLINK2 conversion and merging orchestration script.

    Args:
        bcf_paths (dict[str, str]): Mapping of SG ID to cloud BCF path.
        output_prefix (str): Cloud prefix for the merged PLINK2 files.
        job_name (str): Name for the Hail Batch job.

    Returns:
        Job: A Hail Batch job object.
    """
    b = get_batch()
    j = b.new_job(name=job_name)

    j.image(config_retrieve(['workflow', 'BcfToPlink_image']))
    j.cpu(8)
    j.memory('highmem')
    j.storage('50G')

    # 1. Stage the orchestration script
    script_path = Path(__file__).parent.parent / 'scripts' / 'vcf_to_plink.py'
    vcf_to_plink_script = b.read_input(str(script_path))

    # 2. Stage all BCFs and indices
    staged_bcfs = {}
    for sg_id, cloud_path in bcf_paths.items():
        staged_bcfs[sg_id] = b.read_input_group(
            bcf=cloud_path,
            csi=f'{cloud_path}.csi'
        )

    # 3. Define output resource group
    j.declare_resource_group(
        output_plink={
            'pgen': '{root}.pgen',
            'pvar': '{root}.pvar',
            'psam': '{root}.psam',
        }
    )

    # 4. Construct local manifest and execute script
    # We build the manifest string in python then echo it to a file in the container
    manifest_data = {
        'manifest': {sg_id: str(resource.bcf) for sg_id, resource in staged_bcfs.items()}
    }
    manifest_json = json.dumps(manifest_data)

    j.command(
        f"""
        set -ex

        echo '{manifest_json}' > local_manifest.json

        python3 {vcf_to_plink_script} \\
            --manifest local_manifest.json \\
            --out-prefix {j.output_plink} \\
            --threads 8
        """
    )

    # 5. Write outputs back to cloud
    b.write_output(j.output_plink, output_prefix)

    return j
