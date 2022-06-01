from colour import Color

class element:
    def __init__(this, element_name, element_start, element_end, element_type):
        this.name = element_name
        this.start = element_start
        this.end = element_end
        this.type = element_type
        this.color = Color('black')
        if this.type == "ARG":
            if (this.element_name.find("Metals") != -1) or (this.element_name.find("Biocides") != -1):
                this.color = Color('#27B2D9')
            else:
                this.color = Color('#2A382E')
        elif this.type == "MGE":
            this.color = Color('#EDF511')
        elif this.type == "KEGG":
            this.color = Color('#C7364E')
        
    def getStart(this):
        return this.start

    def getElementInfo(this):
        return (this.name, this.start, this.end, this.color)