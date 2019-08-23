from itcpack import readitc

fixed_diff_requirement=0.1e-6
relative_diff_requirement=0.1

def test_noisy_exothermic_bigheat():
    itcdata=readitc.readitc("datafiles/040419b.itc")
    itcdata.create_titration_dataset()

    #first injection
    trueheat=-2.435e-6
    xsheat_diff=itcdata.xs_heats[0]-trueheat
    assert abs(xsheat_diff)<fixed_diff_requirement \
           or abs(xsheat_diff/trueheat)<relative_diff_requirement

    #second injection
    trueheat=-7.241e-6
    xsheat_diff=itcdata.xs_heats[1]-trueheat
    assert abs(xsheat_diff)<fixed_diff_requirement\
           or abs(xsheat_diff/trueheat)<relative_diff_requirement

def test_exothermic_slowbaselinedecay_injection1_061810b():
    itcdata=readitc.readitc("datafiles/061810b.itc")
    itcdata.create_titration_dataset()
    #first injection
    trueheat=-6.375e-6
    print(itcdata.xs_heats)
    xsheat_diff=itcdata.xs_heats[0]-trueheat
    assert abs(xsheat_diff)<fixed_diff_requirement \
           or abs(xsheat_diff/trueheat)<relative_diff_requirement

def test_exothermic_slowbaselinedecay_injection2_061810b():
    itcdata=readitc.readitc("datafiles/061810b.itc")
    itcdata.create_titration_dataset()
    #second injection
    trueheat=-2.2386e-5
    xsheat_diff=itcdata.xs_heats[1]-trueheat
    assert abs(xsheat_diff)<fixed_diff_requirement\
           or abs(xsheat_diff/trueheat)<relative_diff_requirement

def test_exothermic_slowkinetics_injection9_061810b():
    itcdata=readitc.readitc("datafiles/061810b.itc")
    itcdata.create_titration_dataset()
    #second injection
    trueheat=-2.3416e-5
    xsheat_diff=itcdata.xs_heats[8]-trueheat
    assert abs(xsheat_diff)<fixed_diff_requirement\
           or abs(xsheat_diff/trueheat)<relative_diff_requirement

def test_exothermic_slowkinetics_easyshape_injection10_061810b():
    itcdata=readitc.readitc("datafiles/061810b.itc")
    itcdata.create_titration_dataset()
    #second injection
    trueheat=-1.8685e-5
    xsheat_diff=itcdata.xs_heats[9]-trueheat
    assert abs(xsheat_diff)<fixed_diff_requirement\
           or abs(xsheat_diff/trueheat)<relative_diff_requirement


def test_smallexothermic_slowkinetics_injection13_061810b():
    itcdata=readitc.readitc("datafiles/061810b.itc")
    itcdata.create_titration_dataset()
    #second injection
    trueheat=-7.7813e-6
    xsheat_diff=itcdata.xs_heats[12]-trueheat
    assert abs(xsheat_diff)<fixed_diff_requirement\
           or abs(xsheat_diff/trueheat)<relative_diff_requirement

def test_smallpositive_dilutionheat_injection18_061810b():
    itcdata=readitc.readitc("datafiles/061810b.itc")
    itcdata.create_titration_dataset()
    #second injection
    trueheat=-7.453e-7
    xsheat_diff=itcdata.xs_heats[17]-trueheat
    assert abs(xsheat_diff)<fixed_diff_requirement\
           or abs(xsheat_diff/trueheat)<relative_diff_requirement

def test_smallpositive_dilutionheat_injection21_061810b():
    itcdata=readitc.readitc("datafiles/061810b.itc")
    itcdata.create_titration_dataset()
    #second injection
    trueheat=9.301e-8
    xsheat_diff=itcdata.xs_heats[20]-trueheat
    assert abs(xsheat_diff)<fixed_diff_requirement\
           or abs(xsheat_diff/trueheat)<relative_diff_requirement
