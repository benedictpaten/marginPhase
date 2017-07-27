import logging
import os
from toil_lib import require, UserError
from toil_marginphase.marginphase_pipeline import main


log = logging.getLogger(__name__)



def test_stub(tmpdir):
    require(True, "stub test")
