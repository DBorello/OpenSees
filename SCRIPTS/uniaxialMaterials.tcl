
# ##############################################################
# Define the data structures & procedures for ElasticPP
# ##############################################################

frame .elasticPP

set ElasticPP(materialID) 0
set ElasticPP(E) 0
set ElasticPP(yieldStrain) 0    

# add ElasticPP the materials menu
$m add command -label ElasticPP -command "SetElasticPP"

set count 0
foreach field {materialID E yieldStrain} {
    label .elasticPP.l$field -text $field
    entry .elasticPP.e$field -textvariable ElasticPP($field) -relief sunken 
    grid .elasticPP.l$field -row $count -column 0 -sticky e
    grid .elasticPP.e$field -row $count -column 1 -sticky ew
    incr count
}

button .elasticPP.ok -text OK -command "doneElasticPP"
grid .elasticPP.ok -row 0 -rowspan 3 -column 2 -sticky nsew

proc SetElasticPP { } {
    global toggleFrame
    global .elasticPP
    global matID

    .elasticPP.ematerialID config -state normal
    set ElasticPP(materialID) [expr $matID + 1]
    .elasticPP.ematerialID delete 0 end
    .elasticPP.ematerialID insert 0 $ElasticPP(materialID)
    .elasticPP.ematerialID config -state disabled

    pack forget $toggleFrame
    set toggleFrame .elasticPP
    pack $toggleFrame -side bottom  -fill x
}


proc doneElasticPP { } {
    global matID
    global ElasticPP

    set matID $ElasticPP(materialID)
    uniaxialMaterial ElasticPP $matID $ElasticPP(E) $ElasticPP(yieldStrain)
    eval uniaxialTest $matID
    set matID $ElasticPP(materialID)

    SetValues
    Reset
}



frame .ePP

set EPP(materialID) 0
set EPP(E) 0
set EPP(posY) 0    
set EPP(negY) 0    

# add EPP the materials menu
$m add command -label EPP -command "SetEPP"

set count 0
foreach field {materialID E posY negY} {
    label .ePP.l$field -text $field
    entry .ePP.e$field -textvariable EPP($field) -relief sunken 
    grid .ePP.l$field -row $count -column 0 -sticky e
    grid .ePP.e$field -row $count -column 1 -sticky ew
    incr count
}

button .ePP.ok -text OK -command "doneEPP"
grid .ePP.ok -row 0 -rowspan 3 -column 2 -sticky nsew

proc SetEPP { } {
    global toggleFrame
    global .ePP
    global matID

    .ePP.ematerialID config -state normal
    set EPP(materialID) [expr $matID + 1]
    .ePP.ematerialID delete 0 end
    .ePP.ematerialID insert 0 $EPP(materialID)
    .ePP.ematerialID config -state disabled

    pack forget $toggleFrame
    set toggleFrame .ePP
    pack $toggleFrame -side bottom  -fill x
}


proc doneEPP { } {
    global matID
    global EPP

    set matID $EPP(materialID)
    uniaxialMaterial EPP $matID $EPP(E) $EPP(posY) $EPP(negY)
    eval uniaxialTest $matID
    set matID $EPP(materialID)

    SetValues
    Reset
}




# ##############################################################
# Define the data structures & procedures for Elastic
# ##############################################################

frame .elastic

set Elastic(materialID) 0
set Elastic(E) 0

# add Elastic the materials menu
$m add command -label Elastic -command "SetElastic"

set count 0
foreach field {materialID E} {
    label .elastic.l$field -text $field
    entry .elastic.e$field -textvariable Elastic($field) -relief sunken 
    grid .elastic.l$field -row $count -column 0 -sticky e
    grid .elastic.e$field -row $count -column 1 -sticky ew
    incr count
}

button .elastic.ok -text OK -command "doneElastic"
grid .elastic.ok -row 0 -rowspan 3 -column 2 -sticky nsew

proc SetElastic { } {
    global toggleFrame
    global .elastic
    global matID

    .elastic.ematerialID config -state normal
    set Elastic(materialID) [expr $matID + 1]
    .elastic.ematerialID delete 0 end
    .elastic.ematerialID insert 0 $Elastic(materialID)
    .elastic.ematerialID config -state disabled

    pack forget $toggleFrame
    set toggleFrame .elastic
    pack $toggleFrame -side bottom  -fill x
}


proc doneElastic { } {
    global matID
    global Elastic

    set matID $Elastic(materialID)
    uniaxialMaterial Elastic $matID $Elastic(E)
    eval uniaxialTest $matID
    set matID $Elastic(materialID)

    SetValues
    Reset
}

# ##############################################################
# Define the data structures & procedures for PathIndependent
# ##############################################################

frame .pathInd

set PathInd(materialId) 0
set PathInd(otherMatId) 0

# add PathIndependent to the materials menu
$m add command -label PathIndependent -command "SetPathInd"

set count 0
foreach field {materialId otherMatId} {
    label .pathInd.l$field -text $field
    entry .pathInd.e$field -textvariable PathInd($field) -relief sunken 
    grid .pathInd.l$field -row $count -column 0 -sticky e
    grid .pathInd.e$field -row $count -column 1 -sticky ew
    incr count
}

button .pathInd.ok -text OK -command "donePathInd"
grid .pathInd.ok -row 0 -rowspan 3 -column 2 -sticky nsew

proc SetPathInd { } {
    global toggleFrame
    global .pathInd
    global matID

    .pathInd.ematerialId config -state normal
    set PathInd(materialId) [expr $matID + 1]
    .pathInd.ematerialId delete 0 end
    .pathInd.ematerialId insert 0 $PathInd(materialId)
    .pathInd.ematerialId config -state disabled

    pack forget $toggleFrame
    set toggleFrame .pathInd
    pack $toggleFrame -side bottom  -fill x
}


proc donePathInd { } {
    global matID
    global PathInd

    set matID $PathInd(materialId)
    uniaxialMaterial PathIndependent $matID $PathInd(otherMatId)
    eval uniaxialTest $matID
    set matID $PathInd(materialId)

    SetValues
    Reset
}



# ##############################################################
# Define the data structures & procedures for Hysteretic
# ##############################################################

frame .hysteretic

set Hysteretic(materialId) 0
set Hysteretic(s1p) 0
set Hysteretic(e1p) 0
set Hysteretic(s2p) 0
set Hysteretic(e2p) 0
set Hysteretic(s3p) 0
set Hysteretic(e3p) 0
set Hysteretic(s1n) 0
set Hysteretic(e1n) 0
set Hysteretic(s2n) 0
set Hysteretic(e2n) 0
set Hysteretic(s3n) 0
set Hysteretic(e3n) 0
set Hysteretic(px) 0
set Hysteretic(py) 0
set Hysteretic(d1) 0
set Hysteretic(d2) 0
set Hysteretic(beta) 0

# add Hysteretic to the materials menu
$m add command -label Hysteretic -command "SetHysteretic"

set count 0
foreach field {materialId s1p e1p s2p e2p s3p e3p s1n e1n s2n e2n s3n e3n px py d1 d2 beta} {
    label .hysteretic.l$field -text $field
    entry .hysteretic.e$field -textvariable Hysteretic($field) -relief sunken 
    grid .hysteretic.l$field -row $count -column 0 -sticky e
    grid .hysteretic.e$field -row $count -column 1 -sticky ew
    incr count
}

button .hysteretic.ok -text OK -command "doneHysteretic"
grid .hysteretic.ok -row 0 -rowspan 3 -column 2 -sticky nsew

proc SetHysteretic { } {
    global toggleFrame
    global .hysteretic
    global matID

    .hysteretic.ematerialId config -state normal
    set Hysteretic(materialId) [expr $matID + 1]
    .hysteretic.ematerialId delete 0 end
    .hysteretic.ematerialId insert 0 $Hysteretic(materialId)
    .hysteretic.ematerialId config -state disabled

    pack forget $toggleFrame
    set toggleFrame .hysteretic
    pack $toggleFrame -side bottom  -fill x
}


proc doneHysteretic { } {
    global matID
    global Hysteretic

    set matID $Hysteretic(materialId)
    #uniaxialMaterial Hysteretic $matID $Hysteretic(s1p) $Hysteretic(e1p) $Hysteretic(s2p) $Hysteretic(e2p) $Hysteretic(s3p) $Hysteretic(e3p) $Hysteretic(s1n) $Hysteretic(e1n) $Hysteretic(s2n) $Hysteretic(e2n) $Hysteretic(s3n) $Hysteretic(e3n) $Hysteretic(px) $Hysteretic(py) $Hysteretic(d1) $Hysteretic(d2) $Hysteretic(beta)
    uniaxialMaterial Hysteretic $matID  112 0.012  112 0.020  50 0.025  -110 -0.012 -78 -0.065 -65 -0.080   0.7   0.5   0.01 0.0  0.2

    eval uniaxialTest $matID
    set matID $Hysteretic(materialId)

    SetValues
    Reset
}


# ##############################################################
# Define the data structures & procedures for Hardening
# ##############################################################

frame .hardening

set Hardening(materialId) 0
set Hardening(E) 0
set Hardening(fy) 0
set Hardening(Hiso) 0
set Hardening(Hkin) 0

# add Hardening to the materials menu
$m add command -label Hardening -command "SetHardening"

set count 0
foreach field {materialId E fy Hiso Hkin} {
    label .hardening.l$field -text $field
    entry .hardening.e$field -textvariable Hardening($field) -relief sunken 
    grid .hardening.l$field -row $count -column 0 -sticky e
    grid .hardening.e$field -row $count -column 1 -sticky ew
    incr count
}

button .hardening.ok -text OK -command "doneHardening"
grid .hardening.ok -row 0 -rowspan 3 -column 2 -sticky nsew

proc SetHardening { } {
    global toggleFrame
    global .hardening
    global matID

    .hardening.ematerialId config -state normal
    set Hardening(materialId) [expr $matID + 1]
    .hardening.ematerialId delete 0 end
    .hardening.ematerialId insert 0 $Hardening(materialId)
    .hardening.ematerialId config -state disabled

    pack forget $toggleFrame
    set toggleFrame .hardening
    pack $toggleFrame -side bottom  -fill x
}


proc doneHardening { } {
    global matID
    global Hardening

    set matID $Hardening(materialId)
    uniaxialMaterial Hardening $matID $Hardening(E) $Hardening(fy) $Hardening(Hiso) $Hardening(Hkin)
    eval uniaxialTest $matID
    set matID $Hardening(materialId)

    SetValues
    Reset
}

# ##############################################################
# Define the data structures & procedures for Concrete01
# ##############################################################

frame .concrete01

set Concrete01(materialId) 0
set Concrete01(fc) 0
set Concrete01(ec) 0
set Concrete01(fu) 0
set Concrete01(eu) 0

# add Concrete01 to the materials menu
$m add command -label Concrete01 -command "SetConcrete01"

set count 0
foreach field {materialId fc ec fu eu} {
    label .concrete01.l$field -text $field
    entry .concrete01.e$field -textvariable Concrete01($field) -relief sunken 
    grid .concrete01.l$field -row $count -column 0 -sticky e
    grid .concrete01.e$field -row $count -column 1 -sticky ew
    incr count
}

button .concrete01.ok -text OK -command "doneConcrete01"
grid .concrete01.ok -row 0 -rowspan 3 -column 2 -sticky nsew

proc SetConcrete01 { } {
    global toggleFrame
    global .concrete01
    global matID

    .concrete01.ematerialId config -state normal
    set Concrete01(materialId) [expr $matID + 1]
    .concrete01.ematerialId delete 0 end
    .concrete01.ematerialId insert 0 $Concrete01(materialId)
    .concrete01.ematerialId config -state disabled

    pack forget $toggleFrame
    set toggleFrame .concrete01
    pack $toggleFrame -side bottom  -fill x
}


proc doneConcrete01 { } {
    global matID
    global Concrete01

    set matID $Concrete01(materialId)
    uniaxialMaterial Concrete01 $matID $Concrete01(fc) $Concrete01(ec) $Concrete01(fu) $Concrete01(eu)
    eval uniaxialTest $matID
    set matID $Concrete01(materialId)

    SetValues
    Reset
}



# ##############################################################
# Define the data structures & procedures for Concrete01
# ##############################################################

frame .fedeasconcrete01

set concrete01(materialId) 0
set concrete01(fc) 0
set concrete01(ec) 0
set concrete01(fu) 0
set concrete01(eu) 0

# add concrete01 to the materials menu
$m add command -label concrete01 -command "Setconcrete01"

set count 0
foreach field {materialId fc ec fu eu} {
    label .fedeasconcrete01.l$field -text $field
    entry .fedeasconcrete01.e$field -textvariable concrete01($field) -relief sunken 
    grid .fedeasconcrete01.l$field -row $count -column 0 -sticky e
    grid .fedeasconcrete01.e$field -row $count -column 1 -sticky ew
    incr count
}

button .fedeasconcrete01.ok -text OK -command "doneconcrete01"
grid .fedeasconcrete01.ok -row 0 -rowspan 3 -column 2 -sticky nsew

proc Setconcrete01 { } {
    global toggleFrame
    global .fedeasconcrete01
    global matID

    .fedeasconcrete01.ematerialId config -state normal
    set concrete01(materialId) [expr $matID + 1]
    .fedeasconcrete01.ematerialId delete 0 end
    .fedeasconcrete01.ematerialId insert 0 $concrete01(materialId)
    .fedeasconcrete01.ematerialId config -state disabled

    pack forget $toggleFrame
    set toggleFrame .fedeasconcrete01
    pack $toggleFrame -side bottom  -fill x
}


proc doneconcrete01 { } {
    global matID
    global concrete01

    set matID $concrete01(materialId)
    uniaxialMaterial concrete01 $matID $concrete01(fc) $concrete01(ec) $concrete01(fu) $concrete01(eu)
    eval uniaxialTest $matID
    set matID $concrete01(materialId)

    SetValues
    Reset
}



# ##############################################################
# Define the data structures & procedures for ElasticPPGap
# ##############################################################

frame .elasticPPGap

set elasticPPGap(materialID) 0
set elasticPPGap(E) 0
set elasticPPGap(yieldStrain) 0    
set elasticPPGap(gap) 0    

# add elasticPPGap the materials menu
$m add command -label ElasticPPGap -command "SetElasticPPGap"

set count 0
foreach field {materialID E yieldStrain gap} {
    label .elasticPPGap.l$field -text $field
    entry .elasticPPGap.e$field -textvariable elasticPPGap($field) -relief sunken 
    grid .elasticPPGap.l$field -row $count -column 0 -sticky e
    grid .elasticPPGap.e$field -row $count -column 1 -sticky ew
    incr count
}

button .elasticPPGap.ok -text OK -command "doneElasticPPGap"
grid .elasticPPGap.ok -row 0 -rowspan 3 -column 2 -sticky nsew

proc SetElasticPPGap { } {
    global toggleFrame
    global .elasticPPGap
    global matID

    .elasticPPGap.ematerialID config -state normal
    set elasticPPGap(materialID) [expr $matID + 1]
    .elasticPPGap.ematerialID delete 0 end
    .elasticPPGap.ematerialID insert 0 $elasticPPGap(materialID)
    .elasticPPGap.ematerialID config -state disabled

    pack forget $toggleFrame
    set toggleFrame .elasticPPGap
    pack $toggleFrame -side bottom  -fill x
}


proc doneElasticPPGap { } {
    global matID
    global elasticPPGap

    set matID $elasticPPGap(materialID)
    uniaxialMaterial ElasticPPGap $matID $elasticPPGap(E) $elasticPPGap(yieldStrain) $elasticPPGap(gap)
    eval uniaxialTest $matID
    set matID $elasticPPGap(materialID)

    SetValues
    Reset
}


