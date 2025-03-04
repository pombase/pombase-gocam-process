use std::{collections::{HashMap, HashSet}, fmt::{self, Display}, io::Read};

extern crate serde_json;
#[macro_use] extern crate serde_derive;

use anyhow::Result;

use petgraph::{graph::{NodeIndex, NodeReferences}, visit::{EdgeRef, IntoNodeReferences}, Graph, Undirected};
use petgraph::visit::Bfs;

use pombase_gocam::{FactId, GoCamRawModel, Individual, IndividualId, IndividualType, ModelId, gocam_parse};

pub type GoCamGraph = Graph::<GoCamNode, GoCamEdge>;

pub type GoCamGene = IndividualType;
pub type GoCamChemical = IndividualType;
pub type GoCamModifiedProtein = IndividualType;
pub type GoCamComponent = IndividualType;
pub type GoCamProcess = IndividualType;
pub type GoCamInput = IndividualType;
pub type GoCamOutput = IndividualType;
pub type GoCamGeneIdentifier = String;

pub type GoCamNodeMap = HashMap<IndividualId, GoCamNode>;

#[derive(Clone, Debug)]
pub struct GoCamModel {
    id: String,
    title: String,
    taxon: String,
    graph: GoCamGraph,
}

impl GoCamModel {
    pub fn new(raw_model: GoCamRawModel) -> GoCamModel {
        let graph = make_graph(&raw_model);

        GoCamModel {
            id: raw_model.id().to_owned(),
            title: raw_model.id().to_owned(),
            taxon: raw_model.taxon().to_owned(),
            graph,
        }
    }

    pub fn graph(&self) -> &GoCamGraph {
        &self.graph
    }

    pub fn node_iterator(&self) -> NodeIterator {
        NodeIterator {
            node_refs: self.graph().node_references(),
        }
    }

    pub fn id(&self) -> &str {
        &self.id
    }

    pub fn title(&self) -> &str {
        &self.title
    }

    pub fn taxon(&self) -> &str {
        &self.taxon
    }
}

pub struct NodeIterator<'a> {
    node_refs: NodeReferences<'a, GoCamNode>,
}

impl<'a> Iterator for NodeIterator<'a> {
    type Item = &'a GoCamNode;

    fn next(&mut self) -> Option<Self::Item> {
        self.node_refs.next().map(|(_, node)| node)
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct GoCamComplex {
    pub id: Option<String>,
    pub label: Option<String>,
    pub has_part_genes: Vec<GoCamGeneIdentifier>,
}

impl GoCamComplex {
    pub fn id(&self) -> &str {
        self.id.as_ref().map(|s| s.as_str()).unwrap_or("UNKNOWN_ID")
    }

    pub fn label(&self) -> &str {
        self.label.as_ref().map(|s| s.as_str()).unwrap_or("UNKNOWN_LABEL")
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub enum GoCamEnabledBy {
    Complex(GoCamComplex),
    Gene(GoCamGene),
    Chemical(GoCamChemical),
    ModifiedProtein(GoCamModifiedProtein),
}

impl GoCamEnabledBy {
    pub fn id(&self) -> &str {
        let maybe_id = match self {
            GoCamEnabledBy::Complex(complex) => &complex.id,
            GoCamEnabledBy::Gene(gene) => &gene.id,
            GoCamEnabledBy::Chemical(chemical) => &chemical.id,
            GoCamEnabledBy::ModifiedProtein(modified_protein) => &modified_protein.id,
        };
        maybe_id.as_ref().map(|s| s.as_str()).unwrap_or("UNKNOWN")
    }

    pub fn label(&self) -> &str {
        let maybe_label = match self {
            GoCamEnabledBy::Complex(complex) => &complex.label,
            GoCamEnabledBy::Gene(gene) => &gene.label,
            GoCamEnabledBy::Chemical(chemical) => &chemical.label,
            GoCamEnabledBy::ModifiedProtein(modified_protein) => &modified_protein.label,
        };
        maybe_label.as_ref().map(|s| s.as_str()).unwrap_or("UNKNOWN")
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub enum GoCamNodeType {
    Unknown,
    Chemical,
    Gene(GoCamGene),
    ModifiedProtein(GoCamModifiedProtein),
    Activity(GoCamEnabledBy),
}

#[derive(Serialize, Deserialize, Clone, Debug)]
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

        let (node_type, enabled_by_type, enabled_by_id, enabled_by_label) = match &self.node_type {
            GoCamNodeType::Unknown => ("unknown", "unknown", "unknown", "unknown"),
            GoCamNodeType::Chemical => ("chemical", "", "", ""),
            GoCamNodeType::Gene(_) => ("gene", "", "", ""),
            GoCamNodeType::ModifiedProtein(_) => ("modified_protein", "", "", ""),
            GoCamNodeType::Activity(enabled_by) => match enabled_by {
                GoCamEnabledBy::Chemical(chem) => ("activity", "chemical", chem.id(), chem.label()),
                GoCamEnabledBy::Gene(gene) => ("activity", "gene", gene.id(), gene.label()),
                GoCamEnabledBy::ModifiedProtein(prot) => ("activity", "modified_protein", prot.id(), prot.label()),
                GoCamEnabledBy::Complex(complex) => ("activity", "complex", complex.id(), complex.label()),
            }
        };
        write!(f, "{}\t{}\t{}\t{}\t", node_type, enabled_by_type, enabled_by_id, enabled_by_label)?;

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
            GoCamNodeType::Gene(_) => "gene",
            GoCamNodeType::ModifiedProtein(_) => "modified_protein",
            GoCamNodeType::Activity(activity) => match activity {
                GoCamEnabledBy::Chemical(_) => "enabled_by_chemical",
                GoCamEnabledBy::Gene(_) => "enabled_by_gene",
                GoCamEnabledBy::ModifiedProtein(_) => "enabled_by_modified_protein",
                GoCamEnabledBy::Complex(_) => "enabled_by_complex",
            }
        }
    }

    pub fn description(&self) -> String {
        if let GoCamNodeType::Activity(ref enabler) = self.node_type {
            format!("{} [enabled by] {}", self.label, enabler.label())
        } else {
            self.label.to_owned()
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

#[derive(Serialize, Deserialize, Clone, Debug)]
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
*/
const PROTEIN_CONTAINING_COMPLEX_ID: &str = "GO:0032991";

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
*/

fn individual_is_complex(individual: &Individual) -> bool {
    has_root_term(individual, PROTEIN_CONTAINING_COMPLEX_ID)
}

fn individual_is_gene(individual: &Individual) -> bool {
    if !has_root_term(individual, CHEBI_CHEMICAL_ENTITY_ID) {
        return false;
    }

    if has_root_term(individual, CHEBI_PROTEIN_ID) {
        return false;
    }

    let Some(individual_type) = get_individual_type(individual)
    else {
        return false;
    };

    if let Some(ref id) = individual_type.id {
        if id.starts_with("CHEBI:") {
            return false;
        }
    }

    true
}

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

fn is_connecting_fact(rel_name: &str) -> bool {
    ["causally upstream of, negative effect",
     "causally upstream of, positive effect",
     "provides input for",
     "directly provides input for",
     "removes input for",
     "constitutively upstream of",
     "directly negatively regulates",
     "directly positively regulates",
     "indirectly negatively regulates",
     "indirectly positively regulates",
     "has small molecular activator",
     "has small molecular inhibitor",
     "is small molecule activator of",
     "is small molecule inhibitor of",
     "input of",
     "has output"]
    .iter().any(|s| rel_name == *s)
}

fn make_nodes(model: &GoCamRawModel) -> GoCamNodeMap {
    // genes that are the object of a has_input or has_output relation
    let mut bare_genes_and_modified_proteins = HashSet::new();

    for fact in model.facts() {
        if fact.property_label == "has input" ||
            fact.property_label == "has output" {
                bare_genes_and_modified_proteins.insert(fact.object.clone());
           }
    }

    let mut node_map = HashMap::new();

    for individual in model.individuals() {
        if individual_is_activity(individual) ||
            bare_genes_and_modified_proteins.contains(&individual.id) ||
            individual_is_chemical(individual) &&
            !individual_is_unknown_protein(individual)
        {
            let Some(individual_type) = get_individual_type(individual)
            else {
                continue;
            };
            let detail =
                if individual_is_chemical(individual) {
                    GoCamNodeType::Chemical
                } else {
                    if bare_genes_and_modified_proteins.contains(&individual.id) {
                        if individual_type.id().starts_with("PR:") {
                            GoCamNodeType::ModifiedProtein(individual_type.clone())
                        } else {
                            GoCamNodeType::Gene(individual_type.clone())
                        }
                    } else {
                        GoCamNodeType::Unknown
                    }
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

    let mut complex_map = HashMap::new();

    for individual in model.individuals() {
        if individual_is_complex(individual) {
            let Some(complex_type) = individual.types.get(0)
            else {
                continue;
            };

            let complex = GoCamComplex {
                id: complex_type.id.clone(),
                label: complex_type.label.clone(),
                has_part_genes: vec![],
            };

            complex_map.insert(individual.id.clone(), complex);
        }
    }

    for fact in model.facts() {
        if fact.property_label == "has part" {
            let Some(complex) = complex_map.get_mut(&fact.subject)
            else {
                continue;
            };

            let object_individual = model.fact_object(fact);

            let Some(complex_part_type) = object_individual.types.get(0)
            else {
                continue;
            };

            let complex_gene = complex_part_type.id().to_owned();
//            eprintln!("{}", complex_gene);
            complex.has_part_genes.push(complex_gene);
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
                        let gene_enabler = GoCamEnabledBy::Gene(object_type.clone());
                        subject_node.node_type = GoCamNodeType::Activity(gene_enabler);
                    }
                    else if object_type_id.starts_with("CHEBI:") {
                        let chemical_enabler = GoCamEnabledBy::Chemical(object_type.clone());
                        subject_node.node_type = GoCamNodeType::Activity(chemical_enabler);
                    }
                    else if object_type_id.starts_with("GO:") || object_type_id.starts_with("ComplexPortal:") {
                        let complex = complex_map.get(&fact.object)
                            .expect(&format!("expected complex {}", fact.object))
                            .to_owned();
                        let complex_enabler = GoCamEnabledBy::Complex(complex);
                        subject_node.node_type = GoCamNodeType::Activity(complex_enabler);
                    }
                    else if object_type_id.starts_with("PR:") {
                        let modified_protein_enabler = GoCamEnabledBy::ModifiedProtein(object_type.clone());
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

    let mut connectons = vec![];

    for fact in model.facts() {
        if is_connecting_fact(fact.property_label.as_str()) {
            let connection = (fact.subject.clone(), fact.property_label.clone(),
                              fact.object.clone());
            connectons.push(connection);
        }
    }

    node_map
}

pub fn make_gocam_model(source: &mut dyn Read) -> Result<GoCamModel> {
    let raw_model = gocam_parse(source)?;

    let model = GoCamModel::new(raw_model);

    Ok(model)
}

fn make_graph(model: &GoCamRawModel) -> GoCamGraph {
    let mut graph = GoCamGraph::new();

    let node_map = make_nodes(model);

    let mut id_map = HashMap::new();

    for node in node_map.values() {
        let idx = graph.add_node(node.to_owned());
        id_map.insert(node.individual_gocam_id.clone(), idx);
    }

    for fact in model.facts() {
        let subject_id = &fact.subject;
        let object_id = &fact.object;

        if let (Some(__), Some(_)) =
              (node_map.get(subject_id), node_map.get(object_id))
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

    let graph = &model.graph().clone().into_edge_type::<Undirected>() ;

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
                "gene"|"enabled_by_gene" => {
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

pub fn get_connected_genes(model: &GoCamModel, min_connected_count: usize) -> HashSet<(String, String)> {
    let mut ret = HashSet::new();

    let graph = &model.graph().clone().into_edge_type::<Undirected>();

    let mut seen_idxs = HashSet::new();

    for idx in graph.node_indices() {
        if seen_idxs.contains(&idx) {
            continue;
        }

        seen_idxs.insert(idx);

        let mut connected_genes: HashSet<String> = HashSet::new();

        let mut bfs = Bfs::new(&graph, idx);

        while let Some(nx) = bfs.next(&graph) {
            seen_idxs.insert(nx);

            let gocam_node = graph.node_weight(nx).unwrap();

            match gocam_node.node_type {
                GoCamNodeType::Gene(ref gene) => {
                    connected_genes.insert(gene.id().to_owned());
                },
                GoCamNodeType::Activity(ref enabler) =>
                match enabler {
                    GoCamEnabledBy::Gene(gene) => {
                        connected_genes.insert(gene.id().to_owned());
                    },
                    GoCamEnabledBy::ModifiedProtein(modified_protein) => {
                        connected_genes.insert(modified_protein.id().to_owned());
                    },
                    GoCamEnabledBy::Complex(complex) => {
                        connected_genes.extend(complex.has_part_genes.clone());
                    },
                    _ => (),
                },
                _ => (),
            }
        }

        if connected_genes.len() >= min_connected_count  {
            for gene in connected_genes {
//                println!("{}", gene);
                ret.insert((model.taxon().to_owned(), gene.to_owned()));
            }
        }
    }

    ret
}

type CytoscapeId = String;

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct CytoscapeNodeData {
    pub id: CytoscapeId,
    pub db_id: String,
    pub display_label: String,
    pub label: String,
    pub enabler_label: String,
    #[serde(rename = "type")]
    pub type_string: String,
}

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct CytoscapeEdgeData {
    pub id: CytoscapeId,
    pub label: String,
    pub source: CytoscapeId,
    pub target: CytoscapeId,
    pub weight: i32,
}

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct CytoscapeNode {
    pub data: CytoscapeNodeData,
}

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct CytoscapeEdge {
    pub data: CytoscapeEdgeData
}

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct CytoscapeElements {
    pub nodes: Vec<CytoscapeNode>,
    pub edges: Vec<CytoscapeEdge>,
}

pub fn model_to_cytoscape(model: &GoCamRawModel) -> CytoscapeElements {
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
                    display_label: label.to_owned(),
                    label: label.to_owned(),
                    enabler_label: "".to_owned(),
                }
            })
        }).collect();

    let elements = CytoscapeElements {
        nodes,
        edges,
    };

     elements
}

fn remove_suffix<'a>(s: &'a str, suffix: &str) -> &'a str {
    match s.strip_suffix(suffix) {
        Some(s) => s,
        None => s
    }
}

pub fn model_to_cytoscape_simple(model: &GoCamModel) -> CytoscapeElements {
    let edges: Vec<_> = model.graph().edge_references()
        .map(|edge_ref| {
            let edge = edge_ref.weight();
            let subject_node = model.graph().node_weight(edge_ref.source()).unwrap();
            let object_node = model.graph().node_weight(edge_ref.target()).unwrap();

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

    let nodes: Vec<_> = model.node_iterator()
        .map(|node| {
            let label = remove_suffix(&node.label, " Spom").to_owned();
            let enabler_label = node.enabler_label();
            let enabler_label =
                if enabler_label.len() > 0 {
                    remove_suffix(enabler_label, " Spom").to_owned()
                } else {
                    "".to_owned()
                };
            let display_label = if enabler_label.len() > 0 {
                enabler_label.clone()
            } else {
                label.clone()
            };
            let db_id = node.db_id().to_owned();
            CytoscapeNode {
                data: CytoscapeNodeData {
                    id: node.individual_gocam_id.clone(),
                    db_id,
                    type_string: node.type_string().to_owned(),
                    display_label,
                    label,
                    enabler_label,
                }
            }
        }).collect();

    let elements = CytoscapeElements {
        nodes,
        edges,
    };

    elements
}

pub fn find_holes(model: &GoCamModel) -> Vec<GoCamNode> {
    let node_iter = model.node_iterator();
    node_iter.filter_map(|node| {
        if node.enabler_label() == "protein" {
            Some(node.clone())
        } else {
            None
        }
    }).collect()
}

pub fn find_detached_genes(model: &GoCamRawModel) -> Vec<(String, String, String)> {
    let mut gene_map = HashMap::new();

    for individual in model.individuals() {
        let Some(individual_type) = get_individual_type(individual)
        else {
            continue;
        };

        if individual_is_gene(individual) {
            gene_map.insert(individual.id.clone(), individual_type);
        }
    }

    for fact in model.facts() {
        gene_map.remove(&fact.object);
        gene_map.remove(&fact.subject);
    }

    gene_map.iter().map(|(k, v)| { (k.to_owned(), v.id().to_owned(), v.label().to_owned()) }).collect()
}
