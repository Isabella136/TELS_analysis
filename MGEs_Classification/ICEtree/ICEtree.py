from collections import deque
import threading

already_found = deque()
threads1 = []
ice_to_add_later = []

lock = threading.Lock()

def orderDict(dict_list):
    if len(dict_list) == 1:
        to_return = deque()
        to_return.append(dict_list[0])
        return to_return
    dictA = dict_list[:int(len(dict_list)//2)]
    dictB = dict_list[int(len(dict_list)//2):]
    dictA = orderDict(dictA)
    dictB = orderDict(dictB)
    to_return = deque()
    while (len(dictA) != 0) or (len(dictB) != 0):
        if len(dictA) == 0:
            to_return.append(dictB[0])
            dictB.popleft()
        elif len(dictB) == 0:
            to_return.append(dictA[0])
            dictA.popleft()
        elif len(dictA[0][1]) > len(dictB[0][1]):
            to_return.append(dictA[0])
            dictA.popleft()
        else:
            to_return.append(dictB[0])
            dictB.popleft()
    return to_return

class Head_Node:
    def __init__(this, similarities, q):
        this.similarities = similarities
        this.q = q
        this.ICE_nodes = dict()
    
    def populate(this):
        while (len(this.q) != 0):
            current = this.q.popleft()
            if current in this.ICE_nodes.keys():
                current

class Node:
    def __init__(this, similar_list):
        this.similar_list = similar_list


class ICEtree:
    def __init__(this, similarities):
        this.similarities = dict(orderDict(list(similarities.copy().items())))
        this.head_nodes = []
        this.prep()
        this

    def add_head_node(this, q):
        node = Head_Node(this.similarities, q)
        #Do stuff
        lock.acquire()
        this.head_nodes.append(node)
        lock.release()
    
    def prep(this):
        q = deque()
        for ice, similar_list in this.similarities.items():
            if ice in already_found:            #already parsed through
                continue
            for similar in similar_list:
                if (similar in already_found):  #not parsed through, but should have been
                    ice_to_add_later.append(ice)
                    already_found.append(ice)
                    break
            if ice not in ice_to_add_later:
                this.prep_rec(ice, q)
                if len(q) != 0:
                    threads1.append(threading.Thread(target=this.add_head_node, args=(q.copy(),)))
                    q.clear()
                
    def prep_rec(this, ice, q):
        q.append(ice)
        if ice in already_found:
            q.clear()
            return None
        already_found.append(ice)
        for similar in this.similarities[ice]:
            if similar in q:
                continue
            this.prep_rec(similar, q)
            if len(q) == 0:
                return None
