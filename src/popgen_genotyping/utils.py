"""
suggested location for any utility methods or constants used across multiple stages
"""

from datetime import datetime
from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve
from cpg_flow.workflow import get_workflow

if TYPE_CHECKING:
    from cpg_utils import Path
    from cpg_flow.targets import SequencingGroup, Dataset, Cohort

DATE_STRING: str = datetime.now().strftime('%y-%m')  # noqa: DTZ005


def get_output_prefix(dataset: 'Dataset', stage_name: str, tmp: bool = False) -> 'Path':
    """
    Standardised output prefix for all stages.
    Format: dataset.[tmp_]prefix() / workflow.name / stage_name / version
    """
    version = config_retrieve(['workflow', 'version'], 'v1')
    prefix = dataset.prefix(category='tmp' if tmp else 'main')
    return prefix / get_workflow().name / stage_name / version
