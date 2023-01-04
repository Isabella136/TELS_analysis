from collections import deque
import threading

already_found = deque()
threads1 = []
ice_to_add_later = []

lock = threading.Lock()

class Node:
    def __init__(this, similarities, q):
        this.similarities = similarities
        this.q = q

class ICEtree:
    def __init__(this, similarities):
        this.similarities = similarities
        this.head_nodes = []
        this.prep()

    def add_head_node(this, q):
        node = Node(this.similarities, q)
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
                    # for node in this.head_nodes:
                    #     if node.find(similar):
                    #         node.add(ice)
                    #         break
                    already_found.append(ice)
                    break
            if ice not in ice_to_add_later:
                this.prep_rec(ice, q)
                threads1.append(threading.Thread(target=this.add_head_node, args=(q.copy(),)))
                q.clear()
                
    def prep_rec(this, ice, q):
        q.append(ice)
        for similar in this.similarities[ice]:
            if similar in q:
                continue
            this.prep_rec(similar, q)
