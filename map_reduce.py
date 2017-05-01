from collections import Counter
from collections import defaultdict
import re
import time


def tokenize(message):
    all_words = re.findall(r'\w+', message)
    return set(all_words)


def word_count(documents):
    return Counter(word
                   for document in documents
                   for word in tokenize(document))


def wc_mapper(document):
    for word in tokenize(document):
        yield (word, 1)


def wc_reducer(word, counts):
    yield (word, sum(counts))


def word_count_map_reduce(documents):
    collector = defaultdict(list)

    for document in documents:
        for word, count in wc_mapper(document):
            collector[word].append(count)

    return [output
            for word, counts in collector.items()
            for output in wc_reducer(word, counts)]

start_time = time.time()
print word_count_map_reduce(["data science", "big data", "science fiction"])
print time.time() - start_time

start_time_1 = time.time()
print word_count(["data science", "big data", "science fiction"])
print time.time() - start_time_1