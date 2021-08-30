# screen through titles

rm(list = ls())

d = read.csv('lit_review/unique_entries_abxdur.csv')

## check major papers are present 
# d$Title[grep('treatment for uncomplicated cellulitis', d$Title)] #Comparison of short-course (5 days) and standard (10 days) treatment for uncomplicated cellulitis
# d$Title[grep('Tedizolid', d$Title)] #NOT INCLUDED# Tedizolid Phosphate vs Linezolid for Treatment of Acute Bacterial Skin and Skin Structure Infections: The ESTABLISH-1 Randomized Trial
# d$Title[grep('Tedizolid for 6 days', d$Title)] #Tedizolid for 6 days versus linezolid for 10 days for acute bacterial skin and skin-structure infection
# d$Title[grep('5-day levofloxacin', d$Title)] #Efficacy of 750-mg, 5-day levofloxacin in the treatment of community-acquired pneumonia caused by atypical pathogens
# d$Title[grep('Antibiotic Treatment in Community', d$Title)] #Duration of Antibiotic Treatment in Community-Acquired Pneumonia: A Multicenter Randomized Clinical Trial
# d$Title[grep('8 vs 15 days', d$Title)] #Comparison of 8 vs 15 days of antibiotic therapy for ventilator-associated pneumonia in adults: a randomized trial
# d$Title[grep('empiric antibiotic therapy', d$Title)] #Short-course empiric antibiotic therapy for patients with pulmonary infiltrates in the intensive care unit. A proposed solution for indiscriminate antibiotic prescription
# d$Title[grep('10 days of amoxycillin-clavulanic', d$Title)] #Three vs. 10 days of amoxycillin-clavulanic acid for type 1 acute exacerbations of chronic obstructive pulmonary disease: a randomised
# d$Title[grep('Comparison of ciprofloxacin', d$Title)] #Comparison of ciprofloxacin (7 days) and trimethoprim-sulfamethoxazole (14 days) for acute uncomplicated pyelonephritis pyelonephritis in women
# d$Title[grep('levofloxacin 750 mg once daily', d$Title)] #A trial of levofloxacin 750 mg once daily for 5 days versus ciprofloxacin 400 mg and/or 500 mg twice daily for 10 days in the treatment of acute pyelonephritis
# d$Title[grep('Ciprofloxacin for 7 days versus 14', d$Title)] #Ciprofloxacin for 7 days versus 14 days in women with acute pyelonephritis
# d$Title[grep('Shortened Antimicrobial Treatment', d$Title)] #Shortened Antimicrobial Treatment for Acute Otitis Media in Young Children
# d$Title[grep('Patients with Complicated Intra-Abdominal Infection', d$Title)] #Patients with Complicated Intra-Abdominal Infection Presenting with Sepsis Do Not Require Longer Duration of Antimicrobial Therapy
# d$Title[grep('Seven Versus 14 Days of Antibiotic Therapy', d$Title)] # Seven Versus 14 Days of Antibiotic Therapy for Uncomplicated Gram-negative Bacteremia
# d$Title[grep('Empirical Antimicrobial Therapy in Patients', d$Title)] ##NOT INCLUDED#  Optimisation of Empirical Antimicrobial Therapy in Patients With Haematological Malignancies and Febrile Neutropenia (How Long Study)

###CLEAN 

#####study design - not RCT ,
notRCT.terms = unique(c('rotocol', 'haracteristics', 'post hoc', 'Quality standard', 'survey', 'Design', 'design', 'Myth', 'cohort', 'Factors',
                        'Ventilator-associated pneumonia: Gearing towards shorter-course therapy', 'New treatment alternatives in bacterial vaginosis',
                        'secondary analysis', 'Optimising trial designs', 'Transient bacteremia due to Mycobacterium avium complex in patients with AIDS',
                        'report', 'Comment', 'titer', 'Obesity Is Not Associated', 'Postpneumonic empyema in childhood', 'xperience',
                        'Therapeutic targets in rhinosinusitis: Infection or inflammation?', 'onsiderations', 'ualitative',
                        'Virulence factor expression patterns in Pseudomonas aeruginosa strains from infants with cystic fibrosis',
                        'Prevention of ventilator-associated pneumonia by oral decontamination: Just another SDD study?',
                        'Pulmonary levels of high-mobility group box 1 during mechanical ventilation and ventilator-associated pneumonia',
                        'The efficiency of combination therapy of non-gonococcal urethritis in men', 'bservation', 'Time to First',
                        'Long-term Escherichia coli asymptomatic bacteriuria among women with diabetes mellitus',
                        'Risk Factors', 'factors', 'tympanometric findings', 'Epidemiological', 'association', 'Medical management',
                        'Predicting', 'high risk', 'Management practices', 'overview', 'etrospective',
                        'Recurrent and Relapsing Peritonitis: Causative Organisms and Response to Treatment',
                        'The clinical significance of isolation of two different organisms from the urine of patients with an indwelling catheter',
                        'Can procalcitonin help us in timing of re-intervention in septic patients after multiple trauma or major surgery?',
                        'Role of antibiotics in meconium aspiration syndrome', 'Fall in ICU mortality due to selective decontamination not yet proven',
                        'Spread of Streptococcus pneumoniae and antibiotic-resistant S. pneumoniae from day-care center attendees to their younger siblings',
                        'Procalcitonin and antibiotic use: imperfect, yet effective', 'Asthma exacerbations 5: Assessment and management of severe asthma in adults in hospital',
                        'Molecular characterization of nasal methicillin-resistant staphylococcus aureus isolates showing increasing prevalence of mupirocin resistance and associated multidrug resistance following attempted decolonization',
                        'Management of lactational breast abscesses', 'Management of methicillin-resistant Staphylococcus aureus bacteremia',
                        'go with the flow!', 'Severity of lyme disease with persistent symptoms. Insights from a double-blind placebo-controlled clinical trial',
                        'Empirical antimicrobial therapy of acute dentoalveolar abscess', 'predictors',
                        'Is preoperative distinction between complicated and uncomplicated acute appendicitis feasible without imaging?',
                        'Effect of fluoroquinolone resistance on 14-day levofloxacin triple and triple plus bismuth quadruple therapy',
                        'Antibiotics for treating chronic osteomyelitis in adults', 'Was it really an "ugly" meta-analysis?',
                        'Paradoxical responses after start of antimicrobial treatment in mycobacterium ulcerans infection',
                        'Controlling acute bacterial maxillary sinusitis', 'Macrolide and nonmacrolide resistance with mass azithromycin distribution',
                        'Selective decontamination of the digestive tract: To stimulate or stifle?', 'Update', 'Erratum',
                        'Treating childhood acute lymphoblastic leukemia without cranial irradiation', 'delayed intensification',
                        'Is antimicrobial therapy needed to manage uncomplicated skin and soft-tissue abscesses?',
                        'Does selective digestive decontamination prevent ventilator-associated pneumonia in trauma patients?',
                        'Established and new treatments of the idiopathic inflammatory myopathies: Dermatomyositis and polymyositis',
                        'Mannose-binding lectin deficiency and acute exacerbations of chronic obstructive pulmonary disease',
                        'Intervention strategies for the rational use of anti-microbials', 'Pooled', 'Adherence',
                        'Duration of treatment for nosocomial pneumonia: Choosing the right assessment criteria',
                        'Treatment approaches for infected hemodialysis vascular catheters', 'case control', 'guideline',
                        'Are pediatricians adhering to principles of judicious antibiotic use for upper respiratory tract infections?',
                        '1997-2006', 'Clinical Data', 'emphasis', 'case-control', 'indications for', 'quality improvement',
                        'In defense of evidence: The continuing saga of selective decontamination of the digestive tract',
                        'Demographic, clinical, and treatment parameters influencing the outcome of acute cystitis',
                        'Leprosy type 1 (reversal) reactions and their management', 'ecommend', 'print', 'Hemodialysis',
                        'Vestibular system in infants after systemic therapy with amikacin', 'Determinants',
                        'Antibiotic therapy of osteoarticular infections in the adult', 'etoposide', 
                        'Spread of Streptococcus pneumoniae and antibiotic-resistant S. pneumoniae from day-care center attendees to their younger siblings',
                        'Antibiotics for treating chronic osteomyelitis in adults', 'Cystic fibrosis: Benefits and clinical outcome',
                        'Non-closure of visceral peritoneum at abdominal hysterectomy', 'feasibility study',
                        'Current and future management of serious skin and skin-structure infections', 'Uncomplicated acute bronchitis',
                        'Disease Course', 'Do Polymicrobial Intra-Abdominal Infections Have Worse Outcomes',
                        'multi-state Markov models', 'Predictors', 'presentation'))

notRCT = unlist(lapply(notRCT.terms, function(x){grep(x, d$Title)}))

length(notRCT)

#####participants
notrequireABX.terms = unique(c( 'elicobacter', 'pylori', 'difficile', 'uberculosis', 'recurrent traumatic anterior shoulder instability',
                                'TB', 'tubercular', 'virus', 'viral', 'Is rheumatoid arthritis caused by an infection?',
                                'Viral', 'COVID', 'HIV', 'nfluenza', 'helminth', 'functional dyspepsia',
                                'Is visual outcome compromised when next day review is omitted after phacoemulsification surgery? a randomised control trial',
                                'Does the inactivation of leukocytes in blood transfusions during and following liver transplantation by gamma-irradiation have an impact on rejection and infection rate?',
                                'Does periodontal care improve glycemic control? The Department of Veterans Affairs Dental Diabetes Study',
                                'ataxia', 'myelo', 'psoriasis', 'sclerosis', 'CMV', 'Back Pain in the Elderly',
                                'Interventions for primary vesicoureteric reflux', 'Acute lobar nephronia is associated with a high incidence of renal scarring in childhood urinary tract infections',
                                'Washout policies in long-term indwelling urinary catheterisation in adults.',
                                'Washout policies in long-term indwelling urinary catheterisation in adults.',
                                'A Randomized, Controlled Trial of Panax quinquefolius Extract (CVT-E002) to Reduce Respiratory Infection in Patients With Chronic Lymphocytic Leukemia',
                                'RCT of montelukast as prophylaxis for upper respiratory tract infections in children',
                                'Comparison of 3-Dimensional and Standard Miniplate Fixation in the Management of Mandibular Fractures',
                                'Buprenorphine implants for treatment of opioid dependence: A randomized controlled trial',
                                'Interventions for replacing missing teeth: antibiotics at dental implant placement to prevent complications.',
                                'Routine resite of peripheral intravenous devices every 3 days did not reduce complications compared with clinically indicated resite: A randomised controlled trial',
                                'Use of the PlasmaJet® system in patients undergoing abdominal lipectomy following massive weight loss: A randomized controlled trial',
                                'Relative influence of antibiotic therapy attributes on physician choice in treating acute uncomplicated pyelonephritis',
                                'Randomized Controlled Trial of Enhanced Recovery Program Dedicated to Elderly Patients after Colorectal Surgery',
                                'lymphocyte', 'lupus', 'encephalopathy', 'bowel', 'tumors',
                                'Does periodontal care improve glycemic control? The Department of Veterans Affairs Dental Diabetes Study',
                                'morphology abnormalities', 'laser', 'ooth', 'QT', 'Novel therapies in vasculitis',
                                'cytokine', 'immunoglobulin', 'vein', 'cryptococcal', 'Yeast', 'yeast',
                                'hearing loss', 'eishmaniasis', 'hemorrhoidectomy', 'hypercalc',
                                'Seroma formation after axillary lymphadenectomy with and without the use of drains',
                                'endometritis', 'amebic liver', 'malaria', 'brucellosis', 'polymorphi',
                                'keratitis', 'chronic pelvic pain syndrome', 'toxoplas', 'Candida',
                                'hidradenitis suppurativa', 'Recombinant', 'various clinical',
                                'Trachoma', 'iardia', 'candid', 'leprosy', 'trachoma',
                                'Buruli', 'onchocerciasis', 'Plasmodium falciparum gametocytemia', 
                                'uberculous', 'uchereria bancrofti', 
                                'ityriasis', 'vasculopathy', 'endophthalmitis', 
                                'immune', 'smokers', 'erythema migrans'))

notrequireABX = unlist(lapply(notrequireABX.terms, function(x){grep(x, d$Title)}))

length(notrequireABX)

#####intervention
notAbxdur.terms = unique(c('cortico', 'rednisolone', 'lactoba', 'urokinase', 'Isoprinosine', 'wound coverage', 'disulfiram',
                           'Effect of clarithromycin on esophageal motility', 'mesh', 'Mesh', 'Tranexamic Acid', 'PCNL', 'Heparin','heparin',
                           'Lactoba', 'probiotic', 'Probiotic', 'ranulocyte', 'selenite', 'incision', 'honey', 'T-tube', 'nisin', 'neopterin',
                           '-CSF', 'zinc', 'Zinc', 'adjuvant', 'Topical', 'bupivacaine', 'mesalazine', 'Carboplatin', 'titanium', 'Thalidomide',
                           'Modern application of high voltage stimulation for enhanced healing of venous crural ulceration', 'Onercept', 'Caffeine',
                           'Management of premature rupture of the membranes after 34 weeks', 'ultraso', 'radio', 'Home', 'tenofovir', 'Idebenone', 
                           'external fixator', 'ivermectin','cisplatin', 'Laparoscopy in children with complicated appendicitis', 'Anidulafungin', 'enzastaurin',
                           'topical', 'conjunctivitis', 'solutions', 'oral rehydration solution', 'anti-TNF', 'clamping', 'umeclidinium', 'PNA FISH', 'pleurodesis',
                           'solution in', 'ointment', 'resolution', 'otic solution', 'vaccine', 'lactin-V', 'Day clinic vs. hospital care', 'Mouthwash', 'Complex Approach',
                           'oral solution', 'Itraconazole', 'tobramycin', 'povidone-iodine', 'steroid', 'extended infusion', '4FDC', 'invasive hip surgery',
                           'A randomized trial of ampicillin versus trimethorprim-sulfamethoxazole for 14 days', 'cetirizine',
                           'Randomized trial of amoxicillin for pneumonia in Pakistan', #duplicated
                           'Nebulised hypertonic saline for cystic fibrosis', 'xabepilone', 
                           'Acute cholecystitis: Early versus delayed cholecystectomy, a multicenter randomized trial (ACDC Study, NCT00447304)',
                           'Randomized trial of automated, electronic monitoring to facilitate early detection of sepsis in the intensive care unit',
                           'altrexone', 'Randomized trial of automated, electronic monitoring to facilitate early detection of sepsis in the intensive care unit',
                           'insulin', 'haemorrhoidectomy', 'antioxidants', 'interferon', 'garlic', 'enoxaparin', 'open and closed', 'male circumcision',
                           'Antibiotic treatment of acute otorrhea through tympanostomy tube: randomized double-blind placebo-controlled study with daily follow-up',
                           'antiseptic', 'hlorhex', 'accharomyces boulardii', 'fungal', 'Inhalational versus intravenous', ' delayed graft function', 'aqueous gel',
                           'calcium', 'Cupping', 'chinacea', 'ventricular drain tunneling', 'oropharyngeal aspiration', 'fluorouracil', 'teflon', 'immunotherapy',
                           'Correlates of clinical failure in ventilator-associated pneumonia: insights from a large, randomized trial', 'education', 'Tamsulosin',
                           'globulin', 'Negative Pressure Wound Therapy', 'Trilaciclib', 'Handling', 'Prodefen', 'Saline Drops', 'hyperphosphatemia', 'kin preparation',
                           'Mupirocin', 'spray', 'ventricular repolarization in healthy subjects', 'lumbar microdiscectomy', 'Kushenin', 'approach for acetabular fractures',
                           'Oncologist', 'laparoscopic', 'resuscitation', 'versus 3-trocar needlescopic', 'dosing', 'early versus delayed', 'deltopectoral approaches',
                           'Antibody', 'cysteamine', 'Open Versus Ultrasound Guided', 'Alfacalcidol', 'lancovutide', 'osteotomy', 'hormone', 'statin', 'subantimicrobial',
                           'Fungal', 'aspergillosis', 'mphotericin', 'dynamic', 'grastim', 'zumab', 'Robotic', 'negative pressure', 'regnan', 'risperidone',
                           'vancomycin levels', 'hour', 'Continuous versus intermittent', 'abeprazole', 'staples', 'microwave', 'epirubicin', 'ornidazole',
                           'Continuous versus Intermittent', 'examethasone', 'Steroid', 'ganciclovir', 'Sirolimus', 'tracheostomy versus', 'thoracoscopic debridement',
                           'drainage', 'Drainage', 'impregnate', 'CRP', 'chemotherapy', 'CHOP', 'bilateral total knee arthroplasty', 'wound variables',
                           'sirolimus', 'cell ', 'Fast Breathing', 'antibody', 'anthracycline', 'lenabasum', 'diet', 'tiotropium', 'hand hygiene', 'Phacoemulsification',
                           'diosmectite', 'acupuncture', 'mycophenolate', 'warfarin', 'cyclosporin', 'Uro-Vaxom', 'omeprazole', 'lamotrigine', 'matrix proteins',
                           'stent removal', 'nails', 'Cholecystectomy', 'dressing', 'aunorubicin','Pelargonium', 'proton pump inhibitors', 'coxib', 'XiangBin',
                           'Envelope', 'echnique', 'gastrostomy', 'Piperacillin-Tazobactam vs Meropenem', 'ogurt', 'fixation', 'fenretinide', 'Chang’an I Recip', 
                           'peripherally inserted central catheters vs', 'Lenalidomide', 'Vancomycin versus cefazolin prophylaxis', 'Tocolysis', 'Cysteamine',
                           'Omadacycline', 'prevent bacterial infection', 'olive-fish', 'paracetamol', 'robotic', 'pseudoephedrine', 'Human Neural Stem Cell Transplantation',
                           'Duloxetine', 'Aerosolantibiotic', 'prucalopride', 'tolterodine', 'cream', 'Agomelatine', 'Incisions', 'forceps-guided', 'cellulose',
                           'Pelargonium Extract', 'catheterization', 'interleukin-2', 'tonsillectomy', 'narodustat', 'end expiratory pressure', 'Peppermint',
                           'exercise', 'radiography', 'rehabilitation', 'yolk', 'idotimod', ' Yoghurt', 'Family', 'extract', 'Alzheimer', 'epidural analgesia',
                           'Is visual outcome compromised when next day review is omitted after phacoemulsification surgery? a randomised control trial',
                           'lock solution', 'talactoferrin', 'hand-hygiene', 'fibrinolysis', 'stimulation', 'conduction studies', 'resectoscopy', 'imipenem vs. cefepime',
                           'Do delayed prescriptions reduce antibiotic use in respiratory tract infections? A systematic review', 'Telithromycin', 'ablation',
                           'High-dose versus standard-dose amoxicillin for acute otitis media', 'Percutaneous Nephrolithotomy', 'vincristine', 'adrenaline',
                           'COPP', 'tanercept', 'imab', 'Infectious complications in patients randomized to receive allogeneic bone marrow or peripheral blood transplantation',
                           'Premature rupture of membranes at 34 to 37 weeks', 'cardiopulmonary bypass', 'Teriparatide', 'Anesthesia', 'Laser', 'Diazoxide',
                           'A randomized trial of home vs hospital intravenous antibiotic therapy in adults with infectious diseases', 'epoetin', 'website',
                           'Enteroaggregative Escherichia coli diarrhea in travelers: response to rifaximin therapy', 'berries', 'labour', 'nail', 'thrombin', 'splinting',
                           'Effect of Dermatology Consultation on Outcomes for Patients With Presumed Cellulitis: A Randomized Clinical Trial',
                           'Tissue Plasminogen', 'formula', 'invasive mechanical ventilation', 'umab', 'plasma', 'hand-disinfec', 'Diet', 'Taxane', 'cefditoren pivoxil',
                           'high-frequency chest wall oscillation', 'Kocher-Langenbeck', 'OM-85 BV', 'herb', 'berry', 'Telephone', 'Antiseptic',
                           'pelargonium sidoides', 'oluble fiber', 'nonpharmacological', 'rhubarb', 'nutraceutical', 'sterilization', 'Syphilis', 'axial movement',
                           'itamin', 'lozenges', 'Inulin', 'infection control', 'nystatin', 'everolimus', 'inhibitor', 'regabalin', 'juzen-taiho-to',
                           'surgical evaluation', 'Detection of Methicillin-Resistant', 'ranitidine', 'Adhesive', 'paint', 'Tackers', 'FEC-based therapy', 'salmeterol',
                           'nevirapine', 'C-reactive', 'perioperative oxygen', 'tube feeding', 'Bikini', 'tacrolimus', 'nformation', 'thermometry',
                           'ventilating tube', 'N-methylthiotetrazole', 'Surgical Treatment', 'IBEST', 'sonide', 'pylorus preservation', 'primary total knee replacement',
                           'stents', 'arthroscopic debridement', 'preoperative warming', 'Microbiome Therapeutic', ' comparison of loops for transurethral resection of the prostate',
                           'intravenous lipids', 'cranberr', 'traditional Chinese', 'Bee product', 'miglustat', 'dose-finding', 'Emergency Ureteroscopy',
                           'Jianpi Yiqi', 'acetylcystein', 'ORS', 'sealant', 'accination', 'upplement', 'carbon', 'EURO-LB02', 'surgery or conservative',
                           'percutaneously inserted central venous catheters ', 'STRESS ULCER PROPHYLAXIS', 'fludarabine', 'Device', 'device', 'teriparatide',
                           'Antibiotic treatment of Chlamydia pneumoniae after acute coronary syndrome', 'transobturator tap', 'metadoxine', 'ataluren',
                           'BB-12', 'basiliximab', 'vital signs', 'shave', 'cyclophosphamide', 'Tripterygium wilfordii', 'ivacaftor', 'levodopa-carbidopa',
                           'local', 'implantable', 'Pentoxifylline', 'shaving', 'late intensification therapy', 'mifepristone', 'vasospasm', 'itazoxanide',
                           'suction', 'Local', 'gauze', 'appendectomy', 'octreotide', 'graph', 'taxel', 'eritoran tetrasodium', 'charcoal', 'Incision Drape',
                           'aerosolized drug delivery', 'intravaginal', 'lavage', 'Methenamine', 'Modes of administration', 'switch', 'deoxyribonuclease',
                           'a single oral dose of moxifloxacin (400 mg and 800 mg) on ventricular repolarization in healthy subjects', 'plasminogen',
                           'immunostimulant', 'PEG-gastropexy', 'Synbiotics', 'enhanced versus conventional', 'terazosin', 'Salpingectomy', 'ressing',
                           'rednisone', 'antibod', 'fixation mode', 'intraoral fluoride', 'Prophylactic Endoscopic',
                           'comparing peripherally inserted central venous catheters', 'Inhibitor', 'DRV', 'cortisone', 'saline', 'alfa', 'Vein Grafts',
                           'subsp. lactis', 'coated', 'mannitol', 'lysate'))

notAbxdur = unlist(lapply(notAbxdur.terms, function(x){grep(x, d$Title)}))

length(notAbxdur)

#####outcome
notClinical.terms = unique(c('cost-effectiveness', 'Cost-effectiveness', 'rebleeding', 'conomic', 'Tissue Concentrations',
                             'hospitalization', 'cost of treatment', 'hospital stay', 'functional outcome', 'bactericidal activity',
                             'Duration of fever', 'Renal tolerability', 'quality of life', 'Cost', 'stent', 'weaning',
                             'Pharmacokinetic', 'cost', 'Comparative toxicities', 'Renal Function Changes', 'Satisfaction',
                             'Hospitalization rates', 'cost/effectivity', 'Cardiovascular', 'stimulation', 'colistin versus colistin plus fosfomycin',
                             'economic', 'Microbiological outcomes', 'Hematologic effects', 'gentamicin plus azithromycin and gemifloxacin plus azithromycin',
                             'discharge', 'electrocardiographic changes', 'kinetics'))
notClinical = unlist(lapply(notClinical.terms, function(x){grep(x, d$Title)}))

length(notClinical)

d = d[- unique(c(notRCT, notrequireABX, notAbxdur, notClinical)),]
d = d[order(d$PMID),]
nrow(d)
# keep = unique(c('weeks', 'week', 'days', 'Days', 'Day', 
#                'rocalcitonin', 'uration', 'course', d$Title)))
# d = d[-keep,]

done = read.csv('~/Desktop/done.csv')[,c("relevancy", "reason_exclude", "comments", "meandur_short", "meandur_long", "fu_dur", "PMID", "Title")]
d = d[,c('PMID', 'Title')]
final = merge(d, done, by = 'PMID', all.x = T)
final = final[,c("relevancy", "reason_exclude", "comments", "meandur_short", "meandur_long", "fu_dur", "PMID", "Title.x")]

add = data.frame(relevancy = c(1, 1, NA, NA), 
                 reason_exclude = c('', '', NA, NA), 
                 comments= c('Surveillance sputum, throat, rectum, axilla, and wounds was assessed at the time of ICU admission, once weekly during the stay, at discharge from the ICU, and, if possible, at 7 days after ICU discharge', 
                             NA, NA, NA), 
                 meandur_short = c('', '', NA, NA), 
                 meandur_long = c('', '', NA, NA), 
                 fu_dur= c('', '', NA, NA),  
                 PMID = c(14522530, 11967597, 29153975, 23403680), 
                 Title.x = c("Effects of selective decontamination of digestive tract on mortality and acquisition of resistant bacteria in intensive care: a randomised controlled trial", 
                             "Selective decontamination of subglottic area in mechanically ventilated patients with multiple trauma", 
                             "Optimisation of empirical antimicrobial therapy in patients with haematological malignancies and febrile neutropenia (How Long study): an open-label, randomised, controlled phase 4 trial",
                             "Tedizolid phosphate vs linezolid for treatment of acute bacterial skin and skin structure infections: the ESTABLISH-1 randomized trial"
                 ))

final = rbind(final, add)

write.csv(final, 'lit_review/after_screen_title_withR.csv')
sum(is.na(final$relevancy))

