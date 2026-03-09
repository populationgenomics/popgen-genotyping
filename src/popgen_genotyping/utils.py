"""
suggested location for any utility methods or constants used across multiple stages
"""

from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING

from cpg_utils import to_path, Path as CPGPath
from cpg_utils.config import config_retrieve
from cpg_flow.workflow import get_workflow
from cpg_flow.inputs import get_multicohort

if TYPE_CHECKING:
    from cpg_flow.targets import Dataset, SequencingGroup, Cohort

DATE_STRING: str = datetime.now().strftime('%y-%m')  # noqa: DTZ005


def get_output_prefix(dataset: 'Dataset', stage_name: str, tmp: bool = False) -> 'CPGPath':
    """
    Standardised output prefix for all stages.
    Format: dataset.[tmp_]prefix() / workflow.name / stage_name / version
    """
    version = config_retrieve(['workflow', 'version'], 'v1')
    prefix = dataset.prefix(category='tmp' if tmp else 'main')
    return prefix / get_workflow().name / stage_name / version


def get_sequencing_group_cohort(sequencing_group: 'SequencingGroup') -> 'Cohort':
    """
    Resolve the cohort a sequencing group belongs to by searching the multi-cohort.
    """
    multicohort = get_multicohort()
    for cohort in multicohort.get_cohorts():
        if sequencing_group.id in cohort.get_sequencing_group_ids():
            return cohort

    raise ValueError(f'Sequencing group {sequencing_group.id} not found in any cohort')


def parse_psam(psam_path: str | Path | CPGPath) -> list[str]:
    """
    Extract sequencing group IDs from a PLINK2 .psam file.
    Expects tab-separated values. Sample IDs are in the first column or #IID/IID.
    """
    ids = []
    with to_path(psam_path).open() as f:
        # Find header index for IID
        header = None
        for line in f:
            if not line.strip():
                continue
            if line.startswith('#'):
                header = line.lstrip('#').split()
                iid_idx = header.index('IID') if 'IID' in header else 0
                continue

            parts = line.split()
            if header:
                ids.append(parts[iid_idx])
            else:
                # No header found yet, assume first column
                ids.append(parts[0])

    return ids
