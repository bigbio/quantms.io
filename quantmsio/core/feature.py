import pandas as pd
from quantmsio.core.mzTab import MzTab
from quantmsio.core.sdrf import SDRFHandler
from quantmsio.utils.pride_utils import clean_peptidoform_sequence
from quantmsio.utils.constants import ITRAQ_CHANNEL,TMT_CHANNELS
from quantmsio.core.common import MSSTATS_MAP,MSSTATS_USECOLS
class Feature(MzTab):
    def __init__(self,mzTab_path,sdrf_path,msstats_in_path):
        super(Feature,self).__init__(mzTab_path)
        self._msstats_in = msstats_in_path
        self._sdrf_path = sdrf_path
        self._ms_runs = self.extract_ms_runs()
        self._protein_global_qvalue_map = self.get_protein_map()
        self._modifications = self.get_modifications()
        self._score_names = self.get_score_names()
        self.experiment_type = SDRFHandler(sdrf_path).get_experiment_type_from_sdrf()
        
    def transform_msstats_in(self):
        cols = pd.read_csv(self._msstats_in,nrows=0).columns
        if self.experiment_type == 'LFQ':
            MSSTATS_USECOLS.add('PrecursorCharge')
            MSSTATS_MAP['PrecursorCharge'] = 'precursor_charge'
        else:
            MSSTATS_USECOLS.add('Charge')
            MSSTATS_MAP['Charge'] = 'precursor_charge'
        nocols = MSSTATS_USECOLS - set(cols)
        msstats = pd.read_csv(self._msstats_in,usecols=list(MSSTATS_USECOLS - set(nocols)))
        for col in nocols:
            if col == "Channel":
                msstats.loc[:, col] = "LFQ"
            else:
                msstats.loc[:, col] = None
        msstats["Reference"] = msstats["Reference"].apply(lambda x: x.split(".")[0])
        msstats.loc[:, "sequence"] = msstats["PeptideSequence"].apply(
            lambda x: clean_peptidoform_sequence(x)
        )
        if self.experiment_type != 'LFQ':
            if "TMT" in self.experiment_type:
                msstats["Channel"] = msstats["Channel"].apply(
                    lambda row: TMT_CHANNELS[self.experiment_type][row - 1]
                )
            else:
                msstats["Channel"] = msstats["Channel"].apply(
                    lambda row: ITRAQ_CHANNEL[self.experiment_type][row - 1]
                )
        msstats.rename(columns=MSSTATS_MAP,inplace=True)
        return msstats