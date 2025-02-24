use std::{collections::{HashMap, HashSet}, fmt::{self, Display}};

extern crate serde_json;
#[macro_use] extern crate serde_derive;

use petgraph::{graph::NodeIndex, visit::{EdgeRef, IntoNodeReferences}, Graph, Undirected};
use petgraph::visit::Bfs;

use pombase_gocam::{FactId, GoCamModel, Individual, IndividualId, IndividualType, ModelId};

pub type GoCamGraph = Graph::<GoCamNode, GoCamEdge>;

pub type GoCamComplex = IndividualType;
pub type GoCamGene = IndividualType;
pub type GoCamChemical = IndividualType;
pub type GoCamModifiedProtein = IndividualType;
pub type GoCamComponent = IndividualType;
pub type GoCamProcess = IndividualType;
pub type GoCamInput = IndividualType;
pub type GoCamOutput = IndividualType;

#[derive(Clone, Debug)]
pub enum GoCamActivity {
    Complex(GoCamComplex),
    Gene(GoCamGene),
    Chemical(GoCamChemical),
    ModifiedProtein(GoCamModifiedProtein),
}

impl GoCamActivity {
    pub fn id(&self) -> &str {
        let maybe_id = match self {
            GoCamActivity::Complex(complex) => &complex.id,
            GoCamActivity::Gene(gene) => &gene.id,
            GoCamActivity::Chemical(chemical) => &chemical.id,
            GoCamActivity::ModifiedProtein(modified_protein) => &modified_protein.id,
        };
        maybe_id.as_ref().map(|s| s.as_str()).unwrap_or("UNKNOWN")
    }

    pub fn label(&self) -> &str {
        let maybe_label = match self {
            GoCamActivity::Complex(complex) => &complex.label,
            GoCamActivity::Gene(gene) => &gene.label,
            GoCamActivity::Chemical(chemical) => &chemical.label,
            GoCamActivity::ModifiedProtein(modified_protein) => &modified_protein.label,
        };
        maybe_label.as_ref().map(|s| s.as_str()).unwrap_or("UNKNOWN")
    }
}

#[derive(Clone, Debug)]
pub enum GoCamNodeType {
    Unknown,
    Chemical,
    Activity(GoCamActivity),
}

#[derive(Clone, Debug)]
pub struct GoCamNode {
    pub individual_gocam_id: IndividualId,
    pub id: String,
    pub label: String,
    pub node_type: GoCamNodeType,
    pub has_input: Vec<GoCamInput>,
    pub has_output: Vec<GoCamOutput>,
    pub located_in: Vec<GoCamComponent>,
    pub occurs_in: Vec<GoCamComponent>,
    pub part_of_process: Option<GoCamProcess>,
}

impl Display for GoCamNode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}\t", self.id)?;
        write!(f, "{}\t", self.label)?;

        let (enabled_by_id, enabled_by_label) = match &self.node_type {
            GoCamNodeType::Unknown => ("unknown", "unknown"),
            GoCamNodeType::Chemical => ("unknown_chemical", "unknown_chemical"),
            GoCamNodeType::Activity(activity) => match activity {
                GoCamActivity::Chemical(chem) => (chem.id(), chem.label()),
                GoCamActivity::Gene(gene) => (gene.id(), gene.label()),
                GoCamActivity::ModifiedProtein(prot) => (prot.id(), prot.label()),
                GoCamActivity::Complex(complex) => (complex.id(), complex.label()),
            }
        };
        write!(f, "{}\t{}\t", enabled_by_id, enabled_by_label)?;

        if let Some(ref part_of_process) = self.part_of_process {
            write!(f, "{}\t", part_of_process.label_or_id())?;
        } else {
            write!(f, "\t")?;
        }
        let has_input_string =
            self.has_input.iter().map(|l| l.label_or_id()).collect::<Vec<_>>().join(",");
        if has_input_string.len() > 0 {
            write!(f, "{}\t", has_input_string)?;
        } else {
            write!(f, "\t")?;
        }
        let has_output_string =
            self.has_output.iter().map(|l| l.label_or_id()).collect::<Vec<_>>().join(",");
        if has_output_string.len() > 0 {
            write!(f, "{}\t", has_output_string)?;
        } else {
            write!(f, "\t")?;
        }
        let occurs_in_string =
            self.occurs_in.iter().map(|l| l.label_or_id()).collect::<Vec<_>>().join(",");
        if occurs_in_string.len() > 0 {
            write!(f, "{}\t", occurs_in_string)?
        } else {
            write!(f, "\t")?;
        }
        let located_in_string =
            self.located_in.iter().map(|l| l.label_or_id()).collect::<Vec<_>>().join(",");
        if located_in_string.len() > 0 {
            write!(f, "{}", located_in_string)?;
        }

        Ok(())
    }
}

impl GoCamNode {
    pub fn type_string(&self) -> &str {
        match &self.node_type {
            GoCamNodeType::Unknown => "unknown",
            GoCamNodeType::Chemical => "chemical",
            GoCamNodeType::Activity(activity) => match activity {
                GoCamActivity::Chemical(_) => "enabled_by_chemical",
                GoCamActivity::Gene(_) => "enabled_by_gene",
                GoCamActivity::ModifiedProtein(_) => "enabled_by_modified_protein",
                GoCamActivity::Complex(_) => "enabled_by_complex",
            }
        }
    }

    pub fn enabler_label(&self) -> &str {
        if let GoCamNodeType::Activity(ref enabler) = self.node_type {
            enabler.label()
        } else {
            ""
        }
    }

    pub fn enabler_id(&self) -> &str {
        if let GoCamNodeType::Activity(ref enabler) = self.node_type {
            enabler.id()
        } else {
            ""
        }
    }

    pub fn db_id(&self) -> &str {
        if let GoCamNodeType::Activity(ref enabler) = self.node_type {
            enabler.id()
        } else {
            &self.id
        }
    }
}

#[derive(Clone, Debug)]
pub struct GoCamEdge {
    pub fact_gocam_id: FactId,
    pub id: String,
    pub label: String,
}

impl Display for GoCamEdge {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.label)?;
        Ok(())
    }
}

const MOLECULAR_FUNCTION_ID: &str = "GO:0003674";
/*
const CELLULAR_COMPONENT_ID: &str = "GO:0032991";
const BIOLOGICAL_PROCESS_ID: &str = "GO:0008150";
const PROTEIN_CONTAINING_COMPLEX_ID: &str = "GO:0032991";
*/
const CHEBI_PROTEIN_ID: &str = "CHEBI:36080";
const CHEBI_CHEMICAL_ENTITY_ID: &str = "CHEBI:24431";

fn has_root_term(individual: &Individual, term_id: &str) -> bool {
    for individual_type in &individual.root_types {
        if let Some(ref individual_type_id) = individual_type.id {
            if individual_type_id == term_id {
                return true;
            }
        }
    }

    false
}

fn individual_is_activity(individual: &Individual) -> bool {
    has_root_term(individual, MOLECULAR_FUNCTION_ID)
}

/*
fn individual_is_component(individual: &Individual) -> bool {
    has_root_term(individual, CELLULAR_COMPONENT_ID)
}

fn individual_is_process(individual: &Individual) -> bool {
    has_root_term(individual, BIOLOGICAL_PROCESS_ID)
}

fn individual_is_complex(individual: &Individual) -> bool {
    has_root_term(individual, PROTEIN_CONTAINING_COMPLEX_ID)
}
*/

fn individual_is_chemical(individual: &Individual) -> bool {
    if !has_root_term(individual, CHEBI_CHEMICAL_ENTITY_ID) {
        return false;
    }

    let Some(individual_type) = get_individual_type(individual)
    else {
        return false;
    };

    if let Some(ref id) = individual_type.id {
        if id.starts_with("CHEBI:") {
            return true;
        }
    }

    false
}

fn get_individual_type(individual: &Individual) -> Option<&IndividualType> {
    individual.types.get(0)
}

fn individual_is_unknown_protein(individual: &Individual) -> bool {
    let Some(individual_type) = get_individual_type(individual)
    else {
        return false;
    };

    if let Some(ref individual_type_id) = individual_type.id {
        if individual_type_id == CHEBI_PROTEIN_ID {
            return true;
        }
    }

    false
}

fn is_gene_id(identifier: &str) -> bool {
    ["PomBase:", "FB:", "UniProtKB:", "MGI:", "WB:", "RGD:", "RefSeq:",
     "Xenbase:", "SGD:", "ZFIN:", "RNAcentral:", "EMAPA:"]
        .iter().any(|s| identifier.starts_with(*s))
}

pub fn make_nodes(model: &GoCamModel) -> HashMap<IndividualId, GoCamNode> {
    let mut node_map = HashMap::new();

    for individual in model.individuals() {
        if individual_is_activity(individual) ||
        individual_is_chemical(individual) &&
        !individual_is_unknown_protein(individual) {
            let Some(individual_type) = get_individual_type(individual)
            else {
                continue;
            };
            let detail =
            if individual_is_chemical(individual) {
                GoCamNodeType::Chemical
            } else {
                GoCamNodeType::Unknown
            };
            let gocam_node = GoCamNode {
                individual_gocam_id: individual.id.clone(),
                id: individual_type.id.clone().unwrap_or_else(|| "NO_ID".to_owned()),
                label: individual_type.label.clone().unwrap_or_else(|| "NO_LABEL".to_owned()),
                node_type: detail,
                has_input: vec![],
                has_output: vec![],
                located_in: vec![],
                occurs_in: vec![],
                part_of_process: None,
            };

            node_map.insert(individual.id.clone(), gocam_node);
        }
    }

    for fact in model.facts() {
        let Some(subject_node) = node_map.get_mut(&fact.subject)
        else {
            continue;
        };

        let object_individual = model.fact_object(fact);
        let Some(object_type) = get_individual_type(object_individual)
        else {
            continue;
        };

        match fact.property_label.as_str() {
            "enabled by" => {
                if let Some(ref object_type_id) = object_type.id {
                    if is_gene_id(object_type_id) {
                        let gene_enabler = GoCamActivity::Gene(object_type.clone());
                        subject_node.node_type = GoCamNodeType::Activity(gene_enabler);
                    }
                    else if object_type_id.starts_with("CHEBI:") {
                        let chemical_enabler = GoCamActivity::Chemical(object_type.clone());
                        subject_node.node_type = GoCamNodeType::Activity(chemical_enabler);
                    }
                    else if object_type_id.starts_with("GO:") || object_type_id.starts_with("ComplexPortal:") {
                        let complex_enabler = GoCamActivity::Complex(object_type.clone());
                        subject_node.node_type = GoCamNodeType::Activity(complex_enabler);
                    }
                    else if object_type_id.starts_with("PR:") {
                        let modified_protein_enabler = GoCamActivity::ModifiedProtein(object_type.clone());
                        subject_node.node_type = GoCamNodeType::Activity(modified_protein_enabler);
                    }
                    else  {
                        eprintln!("can't handle enabled by object: {} - {}", object_type_id, object_individual.id);
                    }
                }
            },
            "has input" => {
                subject_node.has_input.push(object_type.clone());
            },
            "has output" => {
                subject_node.has_output.push(object_type.clone());
            },
            "located in" => {
                subject_node.located_in.push(object_type.clone());
            },
            "occurs in" => {
                subject_node.occurs_in.push(object_type.clone());
            },
            "part of" => {
                subject_node.part_of_process = Some(object_type.clone());
            },
            &_ => {
                // eprintln!("ignoring rel from fact: {} {}", fact.property_label, fact.id());
            }
        }
    }

    node_map
}

pub fn make_graph(model: &GoCamModel) -> GoCamGraph {
    let mut graph = GoCamGraph::new();

    let temp_nodes = make_nodes(model);

    let mut id_map = HashMap::new();

    for node in temp_nodes.values() {
        let idx = graph.add_node(node.to_owned());
        id_map.insert(node.individual_gocam_id.clone(), idx);
    }

    for fact in model.facts() {
        let subject_id = &fact.subject;
        let object_id = &fact.object;

        if let (Some(__), Some(_)) =
              (temp_nodes.get(subject_id), temp_nodes.get(object_id))
        {
            let subject_idx = id_map.get(subject_id).unwrap();
            let object_idx = id_map.get(object_id).unwrap();

            let edge = GoCamEdge {
                fact_gocam_id: fact.id(),
                id: fact.property.clone(),
                label: fact.property_label.clone(),
            };

            graph.add_edge(*subject_idx, *object_idx, edge);
        }
    }

    graph
}

pub struct GoCamStats {
    pub id: ModelId,
    pub total_genes: usize,
    pub total_complexes: usize,
    pub max_connected_activities: usize,
    pub total_connected_activities: usize,
    pub number_of_holes: usize,
}

pub fn get_stats(model: &GoCamModel) -> GoCamStats {
    let number_of_holes = find_holes(model).len();

    let graph = &make_graph(&model).into_edge_type::<Undirected>() ;

    let mut seen_idxs = HashSet::new();

    let mut total_genes = 0;
    let mut total_complexes = 0;
    let mut max_connected_activities = 0;
    let mut total_connected_activities = 0;

    for idx in graph.node_indices() {
        if seen_idxs.contains(&idx) {
            continue;
        }

        let mut connected_activities = 0;

        seen_idxs.insert(idx);

        let mut bfs = Bfs::new(&graph, idx);

        let mut inc_counts = |nx: NodeIndex| {
            let gocam_node = graph.node_weight(nx).unwrap();
            let ntype = gocam_node.type_string();
            if let GoCamNodeType::Activity(_) = gocam_node.node_type {
                connected_activities += 1;
            }
            match ntype {
                "enabled_by_gene" => {
                    total_genes += 1;
                },
                "enabled_by_complex" => {
                    total_complexes += 1;
                }
                _ => (),
            }
        };

        while let Some(nx) = bfs.next(&graph) {
            seen_idxs.insert(nx);
            inc_counts(nx);
        }

        if connected_activities > 1 {
            if connected_activities > max_connected_activities {
                max_connected_activities = connected_activities;
            }

            total_connected_activities += connected_activities;
        }

    }

    GoCamStats {
        id: model.id().to_owned(),
        total_genes,
        total_complexes,
        max_connected_activities,
        total_connected_activities,
        number_of_holes,
    }
}

type CytoscapeId = String;

#[derive(Deserialize, Serialize, Debug, Clone)]
struct CytoscapeNodeData {
    id: CytoscapeId,
    db_id: String,
    label: String,
    #[serde(rename = "type")]
    type_string: String,
}

#[derive(Deserialize, Serialize, Debug, Clone)]
struct CytoscapeEdgeData {
    id: CytoscapeId,
    label: String,
    source: CytoscapeId,
    target: CytoscapeId,
    weight: i32,
}

#[derive(Deserialize, Serialize, Debug, Clone)]
struct CytoscapeNode {
    data: CytoscapeNodeData,
}

#[derive(Deserialize, Serialize, Debug, Clone)]
struct CytoscapeEdge {
    data: CytoscapeEdgeData
}

pub fn model_to_cytoscape(model: &GoCamModel) -> String {
    let mut seen_nodes = HashSet::new();

    let edges: Vec<_> = model.facts()
        .map(|fact| {
            seen_nodes.insert(fact.subject.clone());
            seen_nodes.insert(fact.object.clone());

            CytoscapeEdge {
                data: CytoscapeEdgeData {
                    id: fact.id(),
                    label: fact.property_label.clone(),
                    source: fact.subject.clone(),
                    target: fact.object.clone(),
                    weight: 0,
                }
            }
        }).collect();

    let nodes: Vec<_> = model.individuals()
        .filter_map(|individual| {
            if !seen_nodes.contains(&individual.id) {
                return None;
            }

            let Some(individual_type) = individual.types.get(0)
            else {
                return None;
            };

            let individual_type = individual_type.to_owned();

            let (Some(ref label), Some(ref id)) = (individual_type.label, individual_type.id)
            else {
                return None;
            };

            Some(CytoscapeNode {
                data: CytoscapeNodeData {
                    id: individual.id.clone(),
                    db_id: id.clone(),
                    type_string: "node".to_owned(),
                    label: label.to_owned(),
                }
            })
        }).collect();

    let nodes_string = serde_json::to_string(&nodes).unwrap();
    let edges_string = serde_json::to_string(&edges).unwrap();

    format!("nodes: {},\nedges: {}", nodes_string, edges_string)
}

fn remove_suffix<'a>(s: &'a str, suffix: &str) -> &'a str {
    match s.strip_suffix(suffix) {
        Some(s) => s,
        None => s
    }
}

pub fn model_to_cytoscape_simple(graph: &GoCamGraph) -> String {
    let edges: Vec<_> = graph.edge_references()
        .map(|edge_ref| {
            let edge = edge_ref.weight();
            let subject_node = graph.node_weight(edge_ref.source()).unwrap();
            let object_node = graph.node_weight(edge_ref.target()).unwrap();

            let (label, source, target) =
                if edge.label == "has input" {
                    ("input of".into(),
                     object_node.individual_gocam_id.clone(),
                     subject_node.individual_gocam_id.clone())
                } else {
                    (edge.label.to_owned(),
                     subject_node.individual_gocam_id.clone(),
                     object_node.individual_gocam_id.clone())
                };

            CytoscapeEdge {
                data: CytoscapeEdgeData {
                    id: edge.fact_gocam_id.clone(),
                    label,
                    source,
                    target,
                    weight: 0,
                }
            }
        }).collect();

    let nodes: Vec<_> = graph.node_references()
        .map(|(_, node)| {
            let enabler_label = node.enabler_label();
            let label =
                if enabler_label.len() > 0 {
                    remove_suffix(enabler_label, " Spom").to_owned()
                } else {
                    node.label.to_owned()
                };
            let db_id = node.db_id().to_owned();
            Some(CytoscapeNode {
                data: CytoscapeNodeData {
                    id: node.individual_gocam_id.clone(),
                    db_id,
                    type_string: node.type_string().to_owned(),
                    label,
                }
            })
        }).collect();

    let nodes_string = serde_json::to_string(&nodes).unwrap();
    let edges_string = serde_json::to_string(&edges).unwrap();

    format!("nodes: {},\nedges: {}", nodes_string, edges_string)
}

pub fn find_holes(model: &GoCamModel) -> Vec<GoCamNode> {
    let node_map = make_nodes(model);
    node_map.into_values().filter_map(|node| {
        if node.enabler_label() == "protein" {
            Some(node)
        } else {
            None
        }
    }).collect()
}
