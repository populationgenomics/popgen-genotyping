"""
suggested location for any utility methods or constants used across multiple stages
"""

from datetime import datetime
from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve
from cpg_flow.workflow import get_workflow
from cpg_flow.inputs import get_multicohort

if TYPE_CHECKING:
    from cpg_utils import Path
    from cpg_flow.targets import Dataset, SequencingGroup, Cohort

DATE_STRING: str = datetime.now().strftime('%y-%m')  # noqa: DTZ005


def get_output_prefix(dataset: 'Dataset', stage_name: str, tmp: bool = False) -> 'Path':
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
