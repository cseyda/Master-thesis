import matplotlib.pyplot as plt
import pymongo

from numpy.random import normal
#"/home/seydanator/Desktop/SVN/egtdata/tweets.germany.gt"
#  2MB
# 57.554
ger_small = "a56c48c3df72452177dce28efd790ddc"

#"/home/seydanator/Desktop/SVN/egtdata/complete_2008-2011.germany.egt"
#140MB
# 1.254.298
ger = "90bda4cc98289c9d3c231127f8253189"

#usa_small="/home/seydanator/Desktop/SVN/egtdata/complete_2008-2011.usa.egt"
#719MB
# 5.976.723
usa_small = "3f87f7a7aa3d8ba083407ff956f1ada1"

#"/home/seydanator/Desktop/SVN/egtdata/tweets_110716_110815_text.us.egt"
#1200MB
# 14.531.622
usa = "32f974a4edb4eb8f1e2aa366fa1259b1"

connection = pymongo.MongoClient()
stats = connection['tweets']['stats']

#e = stats.find_one({"_id":usa_small})["words_histogram"]
e = stats.find_one({"_id":ger})["bigrams_histogram"]
#gaussian_numbers = normal(size=1000)

plt.hist(e, bins= 100, normed=True)
#plt.title("Gaussian Histogram")
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.show()
