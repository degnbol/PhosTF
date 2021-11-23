#--------------------------------------------------------------------------------
#
#                                 Likelihood module version 1.0
# date: 2013/06/05
#
#--------------------------------------------------------------------------------
import math

class CConvEntry:
    pass

def Count(FCondi, l):
    return sum(map(FCondi, l))

def mean(l):
    sum = 0
    for v in l:
        sum += v
    return float(sum) / len(l)

def ExtendBin(predictions, func_score, lower_bound, upper_bound):
    Bin = []
    for pred in predictions:
        if lower_bound <= func_score(pred) and func_score(pred) <= upper_bound:
            Bin.append(pred)
    return Bin

def IsUpPeak(conv_tbl, i):
    if conv_tbl[i-1].L < conv_tbl[i].L and conv_tbl[i].L >= conv_tbl[i+1].L:
        return True
    else:
        return False
    
def IsDownPeak(conv_tbl, i):
    if conv_tbl[i-1].L >= conv_tbl[i].L and conv_tbl[i].L < conv_tbl[i+1].L:
        return True
    else:
        return False

def Smooth3(conv_tbl, i, num_pos, num_neg, func_score, VAZD):
    num_predictions = num_pos + num_neg
    predictions_bin = list(set(conv_tbl[i].predictions_bin).union(set(conv_tbl[i+1].predictions_bin)).union(set(conv_tbl[i+2].predictions_bin)))
    predictions_bin.sort(key = lambda(x):x.kz)
    predictions_bin.sort(key = func_score, reverse = True)
    num_pos_bin = Count(lambda(x):x.kz == 'k', predictions_bin)
    num_neg_bin = Count(lambda(x):x.kz == 'z', predictions_bin)

    conv_entry = CConvEntry()
    conv_entry.predictions_bin = predictions_bin
    conv_entry.score = mean(map(func_score, predictions_bin))
    conv_entry.score_lower_bound = func_score(predictions_bin[-1])
    conv_entry.score_upper_bound = func_score(predictions_bin[0])

    conv_entry.TPR = float(num_pos_bin) / num_pos
    conv_entry.FPR= float(num_neg_bin) / num_neg
    conv_entry.PPV = float(num_pos_bin) / len(predictions_bin)
    conv_entry.FDR = float(num_neg_bin) / len(predictions_bin)
    conv_entry.L = (float(num_pos_bin) / num_pos) / (float(num_neg_bin + VAZD * (float(len(predictions_bin))/ num_predictions)) / (num_neg + VAZD))
    conv_entry.num_pos = num_pos_bin
    conv_entry.num_neg = num_neg_bin

    conv_tbl.pop(i+2)
    conv_tbl.pop(i+1)
    conv_tbl.pop(i)
    conv_tbl.insert(i, conv_entry)

def MergePoints(conv_tbl, i, num_pos, num_neg, func_score, VAZD):
    num_predictions = num_pos + num_neg
    predictions_bin = list(set(conv_tbl[i].predictions_bin).union(set(conv_tbl[i+1].predictions_bin)))
    predictions_bin.sort(key = lambda(x):x.kz)
    predictions_bin.sort(key = func_score, reverse = True)
    num_pos_bin = Count(lambda(x):x.kz == 'k', predictions_bin)
    num_neg_bin = Count(lambda(x):x.kz == 'z', predictions_bin)

    conv_entry = CConvEntry()
    conv_entry.predictions_bin = predictions_bin
    conv_entry.score = mean(map(func_score, predictions_bin))
    conv_entry.score_lower_bound = func_score(predictions_bin[-1])
    conv_entry.score_upper_bound = func_score(predictions_bin[0])

    conv_entry.TPR = float(num_pos_bin) / num_pos
    conv_entry.FPR= float(num_neg_bin) / num_neg
    conv_entry.PPV = float(num_pos_bin) / len(predictions_bin)
    conv_entry.FDR = float(num_neg_bin) / len(predictions_bin)
    conv_entry.L = (float(num_pos_bin) / num_pos) / (float(num_neg_bin + VAZD * (float(len(predictions_bin))/ num_predictions)) / (num_neg + VAZD))
    conv_entry.num_pos = num_pos_bin
    conv_entry.num_neg = num_neg_bin

    conv_tbl.pop(i)
    conv_tbl.pop(i)
    conv_tbl.insert(i, conv_entry)

    
def SmoothUpPeak(conv_tbl, num_pos, num_neg, func_score, VAZD):
    if len(conv_tbl) < 2:
        #logging.info("Number of entry in conversion table is less than 2")
        return False

    has_up_peak = False
    
    if conv_tbl[-2].L < conv_tbl[-1].L:
        has_up_peak = True
        MergePoints(conv_tbl, len(conv_tbl)-2, num_pos, num_neg, func_score, VAZD)
        i = len(conv_tbl)-5
    else:
        i = len(conv_tbl)-3
        
    if i < 0:
        return has_up_peak

    if len(conv_tbl) >= 3:
        while True:
            if IsUpPeak(conv_tbl, i+1):
                has_up_peak = True
                Smooth3(conv_tbl, i, num_pos, num_neg, func_score, VAZD)
                i -= 3
            else:
                i -= 1
            if i < 0:
                break
        
    return has_up_peak

def SmoothDownPeak(conv_tbl, num_pos, num_neg, func_score, VAZD):
    if len(conv_tbl) < 2:
        #logging.info("Number of entry in conversion table is less than 2")
        return False

    has_down_peak = False

    if len(conv_tbl) >= 3:
        i = len(conv_tbl)-3
        while True:
            if IsDownPeak(conv_tbl, i+1):
                has_down_peak = True
                Smooth3(conv_tbl, i, num_pos, num_neg, func_score, VAZD)
                i -= 3
            else:
                i -= 1
            if i < 0:
                break

    if len(conv_tbl) >= 2:
        if conv_tbl[1].L > conv_tbl[0].L:
            has_down_peak = True
            MergePoints(conv_tbl, 0, num_pos, num_neg, func_score, VAZD)
        
    return has_down_peak

'''
def SmoothReverse(conv_tbl, num_pos, num_neg, func_score, VAZD):
    if len(conv_tbl) < 2:
        logging.info("Number of entry in conversion table is less than 2")
        return False

    has_reverse = False
    
    i = len(conv_tbl)-2
    while i >= 0:
        if conv_tbl[i].L < conv_tbl[i+1].L:
            has_reverse = True
            MergePoints(conv_tbl, i, num_pos, num_neg, func_score, VAZD)
        else:
            i -= 1


    return has_reverse

'''



def RemoveReduncy(conv_tbl, num_pos, num_neg, func_score, VAZD):
    if len(conv_tbl) < 2:
        #logging.info("Number of entry in conversion table is less than 2")
        return False

    has_same_points = False
    
    i = len(conv_tbl)-2
    while True:
        if conv_tbl[i].score == conv_tbl[i+1]:
            has_same_points = True
            MergePoints(conv_tbl, i, num_pos, num_neg, func_score, VAZD)
            i -= 2
        else:
            i -= 1
        if i < 0:
            break

    return has_same_points
                   

def LocalSmooth(conv_tbl, num_pos, num_neg, func_score, VAZD):
    has_same_points = True

    while has_same_points:
        has_same_points = RemoveReduncy(conv_tbl, num_pos, num_neg, func_score, VAZD)
    
    has_up_peak = True
    has_down_peak = True

    conv_tbl.sort(key=lambda(x):x.score, reverse = True)
    
    while True:
        has_up_peak = SmoothUpPeak(conv_tbl, num_pos, num_neg, func_score, VAZD)
        if has_up_peak == False and has_down_peak == False:
            break
        has_down_peak = SmoothDownPeak(conv_tbl, num_pos, num_neg, func_score, VAZD)
        if has_up_peak == False and has_down_peak == False:
            break
'''
        
def LocalSmooth2(conv_tbl, num_pos, num_neg, func_score, VAZD):
    has_same_points = True

    while has_same_points:
        has_same_points = RemoveReduncy(conv_tbl, num_pos, num_neg, func_score, VAZD)
    
    has_reverse = True

    conv_tbl.sort(key=lambda(x):x.score, reverse = True)
    
    while has_reverse:
        has_reverse = SmoothReverse(conv_tbl, num_pos, num_neg, func_score, VAZD)

'''

'''
def WriteConversionTableBin(path_conversion_table, conv_tbl):
    f = open(path_conversion_table, 'w')
    print >> f, "Score\tLower bound\tUpper bound\tLikelihood\tTPR\tFPR\tPPV\tFDR"
    for conv_entry in conv_tbl:
        print >> f, "%.5f\t%.5f\t%.5f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f" % \
              (conv_entry.netphorest_prob, conv_entry.score_lower_bound, conv_entry.score_upper_bound, conv_entry.L, conv_entry.TPR, conv_entry.FPR, conv_entry.PPV, conv_entry.FDR)
    f.close()
'''

def WriteConversionTableBin(path_conversion_table, conv_tbl):
    f = open(path_conversion_table, 'w')
    print >> f, "Score\tLower bound\tUpper bound\tLikelihood\tTPR\tFPR\tPPV\tFDR\tNo. positives\tNo. negatives"
    for conv_entry in conv_tbl:
        print >> f, "%.5f\t%.5f\t%.5f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d" % \
              (conv_entry.score, conv_entry.score_lower_bound, conv_entry.score_upper_bound, conv_entry.L, conv_entry.TPR, conv_entry.FPR, conv_entry.PPV, conv_entry.FDR, conv_entry.num_pos, conv_entry.num_neg)
    f.close()


def ReadConversionTableBin(path_conversion_table):
    conv_tbl = []

    f = open(path_conversion_table)
    for line in f.readlines()[1:]:
        tokens = line.split()
        conv_entry = CConvEntry()
        conv_entry.score = float(tokens[0])
        conv_entry.score_lower_bound = float(tokens[1])
        conv_entry.score_upper_bound = float(tokens[2])
        conv_entry.L = float(tokens[3])
        conv_tbl.append(conv_entry)
        
    f.close()
    
    conv_tbl.sort(key=lambda(x):x.score, reverse=True)
    
    return conv_tbl

    
def WriteConversionTableFDR(path_conversion_table, conv_tbl):
    f = open(path_conversion_table, 'w')
    print >> f, "Score\tFDR"
    for conv_entry in conv_tbl:
        print >> f, "%.5f\t%.3f" % (conv_entry.netphorest_prob, conv_entry.FDR)
    f.close()
    
def GenerateLikelihoodConversionTbl(predictions, num_pos, num_neg, func_score, VAZD):
    bin_size = int(len(predictions)/math.sqrt(num_pos))
    #bin_size = int(float(len(predictions))/num_pos)
    #bin_size = num_pos
    
    conv_tbl = []

    predictions.sort(key = lambda(x):x.kz)
    predictions.sort(key = func_score, reverse = True)
    
    for i in range(0, len(predictions)-bin_size+1):
        #predictions_bin = predictions[i:i+bin_size]
        predictions_bin = ExtendBin(predictions, func_score, func_score(predictions[i+bin_size-1])-0.0001, func_score(predictions[i])+0.0001)
        score = mean(map(func_score, predictions_bin))
        num_pos_bin = Count(lambda(x):x.kz == 'k', predictions_bin)
        num_neg_bin = Count(lambda(x):x.kz == 'z', predictions_bin)
        conv_entry = CConvEntry()
        conv_entry.predictions_bin = predictions_bin
        conv_entry.score = score
        conv_entry.score_lower_bound = func_score(predictions_bin[-1])
        conv_entry.score_upper_bound = func_score(predictions_bin[0])
        conv_entry.TPR = float(num_pos_bin) / num_pos
        conv_entry.FPR= float(num_neg_bin) / num_neg
        conv_entry.PPV = float(num_pos_bin) / len(predictions_bin)
        conv_entry.FDR = float(num_neg_bin) / len(predictions_bin)
        conv_entry.L = (float(num_pos_bin) / num_pos) / (float(num_neg_bin + VAZD * (float(len(predictions_bin))/ len(predictions))) / (num_neg + VAZD))
        conv_entry.num_pos = num_pos_bin
        conv_entry.num_neg = num_neg_bin
        conv_tbl.append(conv_entry)

    return conv_tbl

def WriteConversionTableBin(path_conversion_table, conv_tbl):
    f = open(path_conversion_table, 'w')
    print >> f, "Score\tLower bound\tUpper bound\tLikelihood\tTPR\tFPR\tPPV\tFDR\tNo. positives\tNo. negatives"
    for conv_entry in conv_tbl:
        print >> f, "%.5f\t%.5f\t%.5f\t%.5f\t%.4f\t%.4f\t%.4f\t%.4f\t%d\t%d" % \
              (conv_entry.score, conv_entry.score_lower_bound, conv_entry.score_upper_bound, conv_entry.L, conv_entry.TPR, conv_entry.FPR, conv_entry.PPV, conv_entry.FDR, conv_entry.num_pos, conv_entry.num_neg)
    f.close()


def ReadConversionTableBin(path_conversion_table):
    conv_tbl = []

    f = open(path_conversion_table)
    for line in f.readlines()[1:]:
        tokens = line.split()
        conv_entry = CConvEntry()
        conv_entry.score = float(tokens[0])
        conv_entry.score_lower_bound = float(tokens[1])
        conv_entry.score_upper_bound = float(tokens[2])
        conv_entry.L = float(tokens[3])
        conv_tbl.append(conv_entry)
        
    f.close()
    
    conv_tbl.sort(key=lambda(x):x.score, reverse=True)
    
    return conv_tbl


def ConvertScore2L(score, conv_tbl):
    lower_limit  = 0.0001
    conv_tbl.sort(key=lambda(x):x.score, reverse=True)
    
    # extrapolation
    if score >= conv_tbl[0].score:
        
        if conv_tbl[0].L < lower_limit:
            return lower_limit
        return conv_tbl[0].L
        '''
        L = conv_tbl[-1].L + (conv_tbl[0].L - conv_tbl[-1].L) * (score - conv_tbl[-1].score) / (conv_tbl[0].score - conv_tbl[-1].score)
        if L < conv_tbl[0].L:
            raise "mis calculation"
        if L < lower_limit:
            return lower_limit
        return L
        '''
        
    if score <= conv_tbl[-1].score:
        
        if conv_tbl[-1].L < lower_limit:
            return lower_limit
        return conv_tbl[-1].L
        '''
        L = conv_tbl[-1].L + (conv_tbl[0].L - conv_tbl[-1].L) * (score - conv_tbl[-1].score) / (conv_tbl[0].score - conv_tbl[-1].score)
        if L > conv_tbl[-1].L:
            #1.2510   1.1791  1.1880  0.0243  0.0236  0.0184
            logging.error("%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f" % (L, conv_tbl[0].L, conv_tbl[-1].L, conv_tbl[0].score, conv_tbl[-1].score, score))
            raise "mis calculation"
        if L < lower_limit:
            return 0
        else:
            return L
        '''
        
    
    for i in range(len(conv_tbl)-1):
        if conv_tbl[i].score >= score and score >= conv_tbl[i+1].score:
            if conv_tbl[i].score == conv_tbl[i+1].score:
                if conv_tbl[i].L < lower_limit:
                    return lower_limit
                return conv_tbl[i].L
            else:
                L = conv_tbl[i+1].L + (conv_tbl[i].L - conv_tbl[i+1].L) * (score - conv_tbl[i+1].score) / (conv_tbl[i].score - conv_tbl[i+1].score)
                if L < lower_limit:
                    return lower_limit
                return L

    #logging.error("Can't calculate likelihood")
    return None