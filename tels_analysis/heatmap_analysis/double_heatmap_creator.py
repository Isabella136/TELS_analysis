from tels_analysis.heatmap_analysis.heatmap_creation.indiv_heatmap import IndivHeatmap

class DoubleHeatmapCreator:
    def __init__(self, organism, drug_mech_dict, other_mech_dict):
        self.organism = organism
        self.drug_mech_dict = drug_mech_dict
        self.other_mech_dict = other_mech_dict
        self.ARG_heatmap_list = [
            IndivHeatmap(self.organism, self.drug_mech_dict, "Drugs"),
            IndivHeatmap(self.organism, self.other_mech_dict, "Metals/Biocides")]
        
    def add_to_maps(self, x_axis, filepath, drug_bool_dict, other_bool_dict):
        self.ARG_heatmap_list[0].add_to_map(x_axis, filepath, drug_bool_dict)
        self.ARG_heatmap_list[1].add_to_map(x_axis, filepath, other_bool_dict)

    def make_maps(
            self, drug_bool_dict, other_bool_dict, drug_class_list, other_class_list):
        drug_matrix, drug_x_axis = (
            self.ARG_heatmap_list[0].make_map(drug_bool_dict, drug_class_list))
        other_matrix, other_x_axis = (
            self.ARG_heatmap_list[1].make_map(other_bool_dict, other_class_list))
        return (drug_matrix, drug_x_axis, other_matrix, other_x_axis)