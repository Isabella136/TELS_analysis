
class element:
    def __init__(this, element_name, element_start, element_end, element_type):
        this.name = element_name
        this.start = element_start
        this.end = element_end
        this.type = element_type
        this.color = #000000
        if this.type == "ARG":
            if (this.element_name.find("Metals") != -1) or (this.element_name.find("Biocides") != -1):
                this.color = #27B2D9
            else:
                this.color = #2A382E
        elif this.type == "MGE":
            this.color = #EDF511
        elif this.type == "KEGG":
            this.color = #C7364E
        
    def getStart(this):
        return this.start