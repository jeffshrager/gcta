import os
import gzip
import sys
import numpy 
from sklearn.model_selection import cross_val_predict
from sklearn import linear_model
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt
import numpy as np

def predict(p1, desc):
    result = numpy.matmul(p1,M) 
    print p1, desc, "result=",result
    return result

def print_vec6(name, vec):
    print('%30s: %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f' % (name, 
            lr1.coef_[0], lr1.coef_[1],
            lr1.coef_[2], lr1.coef_[3],
            lr1.coef_[4], lr1.coef_[5]
            ))
def print_vec(name, vec, len):
    print ('%30s:' %name), 
    for i in range(0,len):
        print('%4.1f ' % vec[i]),
    print

if len(sys.argv) != 2:
    print "Please supply filename of input data, must be gzipped"
    sys.exit()

inputfile = sys.argv[1]
filename, file_extension = os.path.splitext(inputfile)
print "filename", filename, "ext", file_extension
if file_extension == ".gz":
    f=gzip.open(inputfile,'rb')
else:
    f = open(inputfile, 'r')
lines = f.readlines()
X = []
X1 = []
X2 = []
X3 = []
X4 = []
X5 = []
X6 = []
X7 = []
X8 = []
X9 = []
X10 = []
X11 = []
Y = []
Y1 = []
Y2 = []
Y3 = []
Y4 = []
Y5 = []
Y6 = []
Y7 = []
Y8 = []
Y9 = []
Y10 = []
Y11 = []
y_pred = []
line_count = 1
drugs = ["BEVACIZUMAB", "TEMOZOLOMIDE", "CABAZITAXEL", "TCAR", "DISATINIB" , "NIVOLUMAB" , "DOXORUBICIN", "DURVALUMAB", "PEMBROLIZUMAB", "VARLILUMAB",  "SORAFENIB" ]

for line in lines:
    if line_count > 1:
        words = line.split()
        bio = int(words[2])
        num = words[3]
        drug = words[4]
        delta_ks = int(words[8])
        #print "biomarker",words[2], "num",words[3], "drug", words[4], "delta-ks",words[7]
        zero17 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        zero7 = [0, 0, 0, 0, 0, 0, 0]

        xnew = zero17[:]

        if (bio & 1) :
            xnew[0] = 1
        if (bio & 2) :
            xnew[1] = 1
        if (bio & 4) :
            xnew[2] = 1
        if (bio & 8) :
            xnew[3] = 1
        if (bio & 16) :
            xnew[4] = 1
        if (bio & 32) :
            xnew[5] = 1

        xdrug = xnew[:6]
        drug_index = drugs.index(drug)
        xnew[drug_index+6] = 1
        #print drug, "index", drug_index

        #for one drug at a time - 6 bit list
        if drug_index == 0:
            X1.append(xdrug)
            Y1.append(delta_ks)
        if drug_index == 1:
            X2.append(xdrug)
            Y2.append(delta_ks)
        if drug_index == 2:
            X3.append(xdrug)
            Y3.append(delta_ks)
        if drug_index == 3:
            X4.append(xdrug)
            Y4.append(delta_ks)
        if drug_index == 4:
            X5.append(xdrug)
            Y5.append(delta_ks)
        if drug_index == 5:
            X6.append(xdrug)
            Y6.append(delta_ks)
        if drug_index == 6:
            X7.append(xdrug)
            Y7.append(delta_ks)
        if drug_index == 7:
            X8.append(xdrug)
            Y8.append(delta_ks)
        if drug_index == 8:
            X9.append(xdrug)
            Y9.append(delta_ks)
        if drug_index == 9:
            X10.append(xdrug)
            Y10.append(delta_ks)
        if drug_index == 10:
            X11.append(xdrug)
            Y11.append(delta_ks)

        #print "delta_ks", delta_ks, output, cab_output, drug
        X.append(xnew)
        Y.append(delta_ks)
        #if line_count < 10:
        #    print "x=",xdrug, "y=",delta_ks

    line_count = line_count + 1

x = np.asarray(X)
x1 = np.asarray(X1)
x2 = np.asarray(X2)
x3 = np.asarray(X3)
x4 = np.asarray(X4)
x5 = np.asarray(X5)
x6 = np.asarray(X6)
x7 = np.asarray(X7)
x8 = np.asarray(X8)
x9 = np.asarray(X9)
x10 = np.asarray(X10)
x11 = np.asarray(X11)
y = np.array(Y)
y1 = np.array(Y1)
y2 = np.array(Y2)
y3 = np.array(Y3)
y4 = np.array(Y4)
y5 = np.asarray(Y5)
y6 = np.asarray(Y6)
y7 = np.asarray(Y7)
y8 = np.asarray(Y8)
y9 = np.asarray(Y9)
y10 = np.asarray(Y10)
y11 = np.asarray(Y11)
lr = linear_model.LinearRegression()
lr1 = linear_model.LinearRegression()
lr2 = linear_model.LinearRegression()
lr3 = linear_model.LinearRegression()
lr4 = linear_model.LinearRegression()
lr5 = linear_model.LinearRegression()
lr6 = linear_model.LinearRegression()
lr7 = linear_model.LinearRegression()
lr8 = linear_model.LinearRegression()
lr9 = linear_model.LinearRegression()
lr10 = linear_model.LinearRegression()
lr11 = linear_model.LinearRegression()

# Train the model using the training sets
lr.fit(x, y)
lr1.fit(x1, y1)
lr2.fit(x2, y2)
lr3.fit(x3, y3)
lr4.fit(x4, y4)
lr5.fit(x5, y5)
lr6.fit(x6, y6)
lr7.fit(x7, y7)
lr8.fit(x8, y8)
lr9.fit(x9, y9)
lr10.fit(x10, y10)
lr11.fit(x11, y11)

predicted = cross_val_predict(lr, x, y, cv=10)
predicted2 = cross_val_predict(lr2, x2, y2, cv=10)
predicted3 = cross_val_predict(lr3, x3, y3, cv=10)

# The coefficients
print_vec('Coefficients all: ', lr.coef_, 17)
print_vec(drugs[0], lr1.coef_, 6)
print_vec(drugs[1], lr2.coef_, 6)
print_vec(drugs[2], lr3.coef_, 6)
print_vec(drugs[3], lr4.coef_, 6)
print_vec(drugs[4], lr5.coef_, 6)
print_vec(drugs[5], lr6.coef_, 6)
print_vec(drugs[6], lr7.coef_, 6)
print_vec(drugs[7], lr8.coef_, 6)
print_vec(drugs[8], lr9.coef_, 6)
print_vec(drugs[9], lr10.coef_, 6)
print_vec(drugs[10], lr11.coef_, 6)

fig, ax = plt.subplots()
ax.scatter(y2, predicted2)
ax.plot([y2.min(), y2.max()], [y2.min(), y2.max()], 'k--', lw=4)
ax.set_xlabel('Measured2')
ax.set_ylabel('Predicted2')
plt.show()

#print "len(y)", len(y), "len(x)", len(x), "len(pred)", len(predicted)
#print "len(y2)", len(y2), "len(x2)", len(x2), "len(pred2)", len(predicted2)
line_count = 0
correct = 0
wrong = 0
for onex in x:
    #y_pred.append(predict(x, "line count=%d y=%s" % (line_count, Y[line_count])))
    try:
        if line_count < 20:
            print onex, y[line_count], predicted[line_count]
        if y[line_count] > 0 and predicted[line_count] > 0:
                correct += 1
        if y[line_count] < 0 and predicted[line_count] < 0:
                correct += 1
        if y[line_count] < 0 and predicted[line_count] > 0:
                wrong += 1
        if y[line_count] > 0 and predicted[line_count] < 0:
                wrong += 1
    except Exception, e:
        pass
    line_count += 1
print "correct", correct, "wrong", wrong, "accuracy %4.1f %%" % (100*float(correct) / float(correct  + wrong))
print "mean rms error", mean_squared_error(y, predicted)
correct = 0
wrong = 0
line_count = 0
for onex in x2:
    #y_pred.append(predict(x, "line count=%d y=%s" % (line_count, Y[line_count])))
    try:
        if line_count < 20:
            print onex, y2[line_count], predicted2[line_count]
        if y2[line_count] > 0 and predicted2[line_count] > 0:
                correct += 1
        if y2[line_count] < 0 and predicted2[line_count] < 0:
                correct += 1
        if y2[line_count] < 0 and predicted2[line_count] > 0:
                wrong += 1
        if y2[line_count] > 0 and predicted2[line_count] < 0:
                wrong += 1
    except Exception, e:
        pass
    line_count += 1
print "correct", correct, "wrong", wrong,"accuracy %6.1f %%" % (100*float(correct) / float(correct  + wrong))
print "mean rms error", mean_squared_error(y, predicted)
for onex in x3:
    try:
        if y3[line_count] > 0 and predicted3[line_count] > 0:
                correct += 1
        if y3[line_count] < 0 and predicted3[line_count] < 0:
                correct += 1
        if y3[line_count] < 0 and predicted3[line_count] > 0:
                wrong += 1
        if y3[line_count] > 0 and predicted3[line_count] < 0:
                wrong += 1
    except Exception, e:
        pass
    line_count += 1
print "correct", correct, "wrong", wrong,"accuracy %6.1f %%" % (100*float(correct) / float(correct  + wrong))
print "mean rms error", mean_squared_error(y, predicted)
