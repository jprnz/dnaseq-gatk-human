def get_sample_dat(record, samples):
    ret = dict()
    gts = record.genotypes
    dps = record.format("AD")

    for i, ind in enumerate(samples):
        ret[ind] = dict()

        dp = sum(dps[i])
        gt = gts[i]
        phase = "|" if gt[2] else "/"

        zyg = None
        if gt[:2] == [-1, -1]:
            pass
        elif gt[:2] == [0,0]:
            zyg = "hom_ref"
        else:
            if len(set(gt[:2])) == 2:
                if 0 not in gt[:2]:
                    zyg = "het_alt"
                else:
                    zyg = "het"
            else:
                zyg = "hom"

        for var in set(gt[:2]):
            ret[ind][var] = dict()
            if var == -1:
                ret[ind][var]["gt"] =  "./."
            else:
                ret[ind][var]["gt"] = phase.join([str(v) for v in gt[:2]])

            ad = int(dps[i][var])
            rd = int(dps[i][0])
            ret[ind][var].update({
                "sample": ind,
                "ad": ad,
                "rd": rd,
                "dp": dp,
                "zyg": zyg})

            if dp > 0:
                ret[ind][var].update({
                    "ap": 100 * ad / dp,
                    "rp": 100 * rd / dp})
            else:
                ret[ind][var].update({ "ap": 0, "rp": 0})

    return ret

