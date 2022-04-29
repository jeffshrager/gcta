import nltk
import csv

# Turns an array of words into training '{contents{word}':True,...} for the  nltk learner
def contentsize(words):
    # UUU This can probably be done in a comprehension
    features = {}
    for word in words:
        features['contains({})'.format(word)] = True
    return features

global train_set, test_set, featuresets, treatments

def init():
    global treatments
    # Import and parse genetic keys
    with open("3705680844-genetics.tsv") as i:
        reader = csv.reader(i, delimiter="\t")
        genetics = {k: v.split() for (k,v) in list(reader)}
    # Import and parse treatment experiences
    with open("3705680844-mrhlog.xls") as i: 
        reader = csv.reader(i, delimiter="\t")
        treatments = [[contentsize(genetics[mkey]),tx,score] for [ingore0,mkey,ignore2,tx,i4,i5,i6,score] in list(reader)]

# Create complete learning set, and both the training and test sets
def dolearn(threshold=10):
    global train_set, test_set, featuresets, treatments
    featuresets = [(elt[0],elt[1]) for elt in treatments if (int(elt[2])>=threshold)]
    print len(featuresets)
    size = int(len(featuresets) * 0.1)
    train_set, test_set = featuresets[size:], featuresets[:size]
    classifier = nltk.NaiveBayesClassifier.train(train_set)
    return (nltk.classify.accuracy(classifier, test_set))

init()
for i in range(-20,20):
    print(str(i)+":"+str(dolearn(i)))

