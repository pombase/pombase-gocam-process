use std::{collections::{HashMap, HashSet}, fmt::{self, Display}};

extern crate serde_json;
#[macro_use] extern crate serde_derive;

use petgraph::{dot::Config, graph::NodeIndex, visit::{EdgeRef, IntoNodeReferences, NodeRef}, Graph, Undirected};
use petgraph::visit::Bfs;

use pombase_gocam::{FactId, GoCamModel, Individual, IndividualId, IndividualType};

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
enum GoCamActivity {
    Complex(GoCamComplex),
    Gene(GoCamGene),
    Chemical(GoCamChemical),
    ModifiedProtein(GoCamModifiedProtein),
}

impl GoCamActivity {
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
        write!(f, "{}\t", self.enabler_label())?;
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
            write!(f, "{}", occurs_in_string)?
        }
        let located_in_string =
            self.located_in.iter().map(|l| l.label_or_id()).collect::<Vec<_>>().join(",");
        if located_in_string.len() > 0 {
            write!(f, "{}\t", located_in_string)?;
        } else {
            write!(f, "\t")?;
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

fn make_nodes(model: &GoCamModel) -> HashMap<IndividualId, GoCamNode> {
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

fn make_graph(model: &GoCamModel) -> GoCamGraph {
/*
    let model_id = model.id();
    let model_title = model.title();
    let model_taxon = model.taxon();
*/
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

            /*
            println!("{}: {} ({}) <- {} -> {} ({})",
                     model.id(),
                     subject_node.label, subject_node.id,
                     fact.property_label,
                     object_node.label, object_node.id);
            */
        }
    }
    use petgraph::dot::Dot;

let dag_graphviz = Dot::with_attr_getters(
    &graph,
    &[Config::NodeNoLabel, Config::EdgeNoLabel],
    &|_, edge| format!("label = \"{}\"", edge.weight().label),
    &|_, node| {
        let enabler_label = node.weight().enabler_label();
        if enabler_label.len() > 0 {
            format!("label = \"{}\"", enabler_label)
        } else {
            format!("label = \"{}\"", node.weight().label)
        }
    },
);

//    println!("{}", Dot::new(&graph));
//    println!("{}", dag_graphviz);

    graph
}

pub fn print_tuples(model: &GoCamModel) {
    let empty = &"".to_owned();
    for fact in model.facts() {
        let subject = model.fact_subject(fact);
        let object = model.fact_object(fact);
        let Some(subject_type) = subject.types.get(0)
        aelse {
            continue;
        };
        let Some(object_type) = object.types.get(0)
        else {
            continue;
        };
        println!("{}\t{}\t{}\t{}\t{}\t{}\t{}",
                 model.id(),
                 model.title(),
                 subject_type.label.as_ref().unwrap_or(empty),
                 subject_type.id.as_ref().unwrap_or(&subject_type.type_string),
                 fact.property_label,
                 object_type.label.as_ref().unwrap_or(empty),
                 object_type.id.as_ref().unwrap_or(&object_type.type_string));
    }
}

pub fn print_stats(model: &GoCamModel) {
    let graph = &make_graph(&model).into_edge_type::<Undirected>() ;

    let mut seen_idxs = HashSet::new();

    let mut total_genes = 0;
    let mut max_connected_genes = 0;
    let mut total_connected_genes = 0;

    for idx in graph.node_indices() {
        if seen_idxs.contains(&idx) {
            continue;
        }

        let mut connected_genes = 0;

        seen_idxs.insert(idx);

        let mut bfs = Bfs::new(&graph, idx);

        let mut inc_counts = |nx: NodeIndex| {
            let gocam_node = graph.node_weight(nx).unwrap();
            let ntype = gocam_node.type_string();
            if ntype.starts_with("enabled_by") {
                total_genes += 1;
                connected_genes += 1;
            }
        };

        while let Some(nx) = bfs.next(&graph) {
            seen_idxs.insert(nx);
            inc_counts(nx);
        }

        if connected_genes > 1 {
            if connected_genes > max_connected_genes {
                max_connected_genes = connected_genes;
            }

            total_connected_genes += connected_genes;
        }
    }

    println!("{}\t{}\t{}\t{}\t{}", model.id(), model.taxon(),
             total_genes, max_connected_genes, total_connected_genes);
}

type CytoscapeId = String;

#[derive(Deserialize, Serialize, Debug, Clone)]
struct CytoscapeNodeData {
    id: CytoscapeId,
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
}

#[derive(Deserialize, Serialize, Debug, Clone)]
struct CytoscapeNode {
    data: CytoscapeNodeData,
}

#[derive(Deserialize, Serialize, Debug, Clone)]
struct CytoscapeEdge {
    data: CytoscapeEdgeData
}

fn model_to_cytoscape(model: &GoCamModel) -> String {
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
        let label = format!("{} ({})", label, id);
        Some(CytoscapeNode {
            data: CytoscapeNodeData {
                id: individual.id.clone(),
                type_string: "node".to_owned(),
                label,
            }
        })
     }).collect();

     let nodes_string = serde_json::to_string(&nodes).unwrap();
     let edges_string = serde_json::to_string(&edges).unwrap();

     format!("nodes: {},\nedges: {}", nodes_string, edges_string)
}

fn model_to_cytoscape_simple(graph: &GoCamGraph) -> String {
     let edges: Vec<_> = graph.edge_references()
         .map(|edge_ref| {
             let edge = edge_ref.weight();
             let subject_node = graph.node_weight(edge_ref.source()).unwrap();
             let object_node = graph.node_weight(edge_ref.target()).unwrap();

             CytoscapeEdge {
                 data: CytoscapeEdgeData {
                     id: edge.fact_gocam_id.clone(),
                     label: edge.label.clone(),
                     source: subject_node.individual_gocam_id.clone(),
                     target: object_node.individual_gocam_id.clone(),
                 }
             }
         }).collect();

     let nodes: Vec<_> = graph.node_references()
         .map(|(_, node)| {
             let enabler_label = node.enabler_label();
             let label =
                 if enabler_label.len() > 0 {
                     enabler_label.to_owned()
                 } else {
                     format!("{} ({})", node.label, node.id)
                 };
             Some(CytoscapeNode {
                 data: CytoscapeNodeData {
                     id: node.individual_gocam_id.clone(),
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
    })
    .collect()
}
