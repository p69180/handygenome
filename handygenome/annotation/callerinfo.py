from handygenome.annotation.annotitem import AnnotItemInfoSingle


class CallerInfo(AnnotItemInfoSingle):
    meta = {
        "ID": "caller_info",
        "Number": "1",
        "Type": "String",
        "Description": "Information about which variant calling program called this variant.",
    }

    @classmethod
    def from_vr(cls, vr):
        return cls.from_vr_base(vr)

    def write(self, vr):
        self.write_base(vr)

    @classmethod
    def from_caller_names(cls, caller_names):
        result = cls.init_nonmissing()
        result['callers'] = tuple(caller_names)
        return result
    

def get_merged_callerinfo(vr_list):
    raw_dict = {'callers': set()}
    for vr in vr_list:
        if CallerInfo.get_annotkey() in vr.info.keys():
            unit_callerinfo = CallerInfo.from_vr(vr)
            raw_dict['callers'].update(unit_callerinfo['callers'])
    raw_dict['callers'] = tuple(raw_dict['callers'])

    result = CallerInfo.init_nonmissing()
    result.update(raw_dict)
    return result
