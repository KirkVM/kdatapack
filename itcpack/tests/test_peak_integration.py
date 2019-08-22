from itcpack import readitc

fixed_diff_requirement=0.1e-6
relative_diff_requirement=0.1

def test_noisy_exothermic_bigheat_data():
    itcdata=readitc.readitc("datafiles/040419b.itc")
    itcdata.create_titration_dataset()

    #first injection
    knownval=-2.435e-6
    xsheat_diff=itcdata.xs_heats[0]-knownval
    assert abs(xsheat_diff)<fixed_diff_requirement \
           or abs(xsheat_diff/knownval)<relative_diff_requirement

    #second injection
    knownval=-7.241e-6
    xsheat_diff=itcdata.xs_heats[1]-knownval
    assert abs(xsheat_diff)<fixed_diff_requirement\
           or abs(xsheat_diff/knownval)<relative_diff_requirement

def test_exothermic_smallheat_data():
    itcdata=readitc.readitc("datafiles/061810b.itc")
    itcdata.create_titration_dataset()

    #first injection
    knownval=-6.375e-6
    print(itcdata.xs_heats)
    xsheat_diff=itcdata.xs_heats[0]-knownval
    assert abs(xsheat_diff)<fixed_diff_requirement \
           or abs(xsheat_diff/knownval)<relative_diff_requirement

    #second injection
    knownval=-2.239e-5
    xsheat_diff=itcdata.xs_heats[1]-knownval
    assert abs(xsheat_diff)<fixed_diff_requirement\
           or abs(xsheat_diff/knownval)<relative_diff_requirement