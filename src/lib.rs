use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

extern crate serde_json;
#[macro_use] extern crate serde_derive;

#[macro_use] extern crate lazy_static;

use itertools::Itertools;

use petgraph::{graph::NodeIndex, Undirected};
use petgraph::visit::Bfs;
use petgraph::visit::EdgeRef;

use pombase_gocam::{GoCamDirection, GoCamGeneIdentifier};
use pombase_gocam::{GoCamComponent, GoCamEnabledBy, GoCamModel,
                    GoCamNode, GoCamNodeOverlap, GoCamProcess,
                    GoCamNodeType, raw::GoCamRawModel, ModelId};
use regex::Regex;

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

pub type GoCamConnectedByCount = HashMap<usize, HashSet<String>>;

pub fn get_connected_genes(model: &GoCamModel)
 -> GoCamConnectedByCount
{
    let mut ret = GoCamConnectedByCount::new();

    let graph = &model.graph().clone().into_edge_type::<Undirected>();

    let mut seen_idxs = HashSet::new();

    for idx in graph.node_indices() {
        if seen_idxs.contains(&idx) {
            continue;
        }

        seen_idxs.insert(idx);

        let mut connected_genes = HashSet::new();
        let mut seen_activities = HashSet::new();

        let mut bfs = Bfs::new(&graph, idx);

        while let Some(nx) = bfs.next(&graph) {
            seen_idxs.insert(nx);

            let gocam_node = graph.node_weight(nx).unwrap();

            match gocam_node.node_type {
                GoCamNodeType::Gene(ref gene) => {
                    connected_genes.insert(gene.id().to_owned());
                },
                GoCamNodeType::Activity(ref enabler) => {
                    seen_activities.insert(gocam_node.individual_gocam_id.clone());
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
                    }
                },
                _ => (),
            }
        }

        for count in 0..=seen_activities.len() {
            ret.entry(count)
               .or_insert_with(HashSet::new)
               .extend(connected_genes.iter().cloned());
        }
    }

    ret
}

type CytoscapeId = String;

#[derive(Deserialize, Serialize, Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct CytoscapeNodeData {
    pub id: CytoscapeId,
    pub display_label: String,
    pub node_id: String,
    pub label: String,
    pub enabler_label: String,
    pub enabler_id: String,
    #[serde(skip_serializing_if="Option::is_none")]
    pub located_in: Option<GoCamComponent>,
    #[serde(skip_serializing_if="BTreeSet::is_empty")]
    pub occurs_in: BTreeSet<GoCamComponent>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub part_of_process: Option<GoCamProcess>,
    pub has_part_genes: BTreeSet<GoCamGeneIdentifier>,
    // this node is in more than one model
    pub is_connecting_node: bool,
    #[serde(rename = "type")]
    pub type_string: String,
    pub model_ids: BTreeSet<ModelId>,
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

#[derive(Serialize, Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct CytoscapeModel {
    pub id: CytoscapeId,
    pub title: String,
}

#[derive(Serialize, Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct CytoscapeModelConnection {
    pub id: CytoscapeId,
    pub label: String,
    pub source: ModelId,
    pub target: ModelId,
    pub enabler_id: String,
    pub has_direction: bool,
}

#[derive(Serialize, Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct CytoscapeModelElements {
    pub nodes: BTreeSet<CytoscapeModel>,
    pub edges: BTreeSet<CytoscapeModelConnection>,
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

            let mut model_ids = BTreeSet::new();
            model_ids.insert(model.id().to_owned());

            Some(CytoscapeNode {
                data: CytoscapeNodeData {
                    id: individual.id.clone(),
                    type_string: "node".to_owned(),
                    display_label: label.to_owned(),
                    label: label.to_owned(),
                    node_id: id.clone(),
                    enabler_label: "".to_owned(),
                    enabler_id: "".to_owned(),
                    located_in: None,
                    occurs_in: BTreeSet::new(),
                    part_of_process: None,
                    has_part_genes: BTreeSet::new(),
                    is_connecting_node: false,
                    model_ids,
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

pub fn model_to_cytoscape_simple(model: &GoCamModel, overlaps: &Vec<GoCamNodeOverlap>) -> CytoscapeElements {
    let mut overlap_set = BTreeSet::new();

    for overlap in overlaps {
        overlap_set.extend(overlap.overlapping_individual_ids.iter().cloned());
    }

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
        .map(|(_, node)| {
            let label = remove_suffix(remove_suffix(&node.label, " Spom"),
                                      " Dmel")
                .to_owned();
            let enabler_label = node.enabler_label();
            let enabler_label =
                if node.enabler_id() == "CHEBI:36080" {
                    "unknown protein".to_owned()
                } else {
                    if enabler_label.len() > 0 {
                        remove_suffix(remove_suffix(enabler_label, " Spom"), " Dmel").to_owned()
                    } else {
                        "".to_owned()
                    }
                };
            let enabler_id = node.enabler_id().to_owned();
            let display_label =
                if enabler_label.len() > 0 {
                    enabler_label.clone()
                } else {
                    label.clone()
                };

            let is_connecting_node =
                 overlap_set.intersection(&node.source_ids).next().is_some();
            let has_part_genes =
                if let GoCamNodeType::Activity(ref enabled_by) = node.node_type {
                    if let GoCamEnabledBy::Complex(complex) = enabled_by {
                        complex.has_part_genes.iter()
                            .map(|id| {
                                if let Some(new_id) = id.split(':').nth(1) {
                                    new_id.to_owned()
                                } else {
                                    id.to_owned()
                                }
                            })
                            .collect()
                    } else {
                        BTreeSet::new()
                    }
                } else {
                    BTreeSet::new()
                };
            let model_ids = node.models.iter().map(|(model_id, _)| model_id.to_owned()).collect();
            CytoscapeNode {
                data: CytoscapeNodeData {
                    id: node.individual_gocam_id.clone(),
                    type_string: node.type_string().to_owned(),
                    display_label,
                    label,
                    node_id: node.node_id.clone(),
                    enabler_label,
                    enabler_id,
                    located_in: node.located_in.clone(),
                    occurs_in: node.occurs_in.clone(),
                    part_of_process: node.part_of_process.clone(),
                    has_part_genes,
                    is_connecting_node,
                    model_ids,
                }
            }
        }).collect();

    let elements = CytoscapeElements {
        nodes,
        edges,
    };

    elements
}

pub fn model_connections_to_cytoscope(overlaps: &Vec<GoCamNodeOverlap>, model_ids_and_titles: &[(String, String)])
   -> CytoscapeModelElements
{
    let mut edges = BTreeSet::new();

    let mut models_in_overlaps = HashSet::new();

    for overlap in overlaps {
        let GoCamNodeType::Activity(ref enabler) = overlap.node_type
        else {
            continue;
        };

        let iter = overlap.models.iter()
            .cartesian_product(overlap.models.iter());

        for (first, second) in iter.into_iter() {
            if first >= second {
                continue;
            }

            models_in_overlaps.insert(first.to_owned());
            models_in_overlaps.insert(second.to_owned());

            let (first_model_id, _, first_direction) = first.to_owned();
            let (second_model_id, _, second_direction) = second.to_owned();

            let (source, target) =
                if first_direction == second_direction {
                    continue;
                } else {
                    if first_direction == GoCamDirection::Incoming ||
                        second_direction == GoCamDirection::Outgoing
                    {
                        (first_model_id, second_model_id).to_owned()
                    } else {
                        (second_model_id, first_model_id).to_owned()
                    }
                };

            let enabler_id = enabler.id().to_owned();

            let edge = CytoscapeModelConnection {
                id: format!("{}-{}", source, target),
                label: overlap.node_label.to_owned(),
                source,
                target,
                enabler_id,
                has_direction: true,
            };
            edges.insert(edge);
        }
    }

    let mut nodes: BTreeSet<_> = models_in_overlaps.iter()
        .map(|(model_id, model_title, _)| {
            CytoscapeModel {
                id: model_id.to_owned(),
                title: model_title.to_owned(),
            }
        })
        .collect();

    for (id, title) in model_ids_and_titles {
        let cytoscape_model = CytoscapeModel {
            id: id.clone(),
            title: title.clone(),
        };

        nodes.insert(cytoscape_model);
    }

    CytoscapeModelElements {
        nodes,
        edges,
    }
}

pub fn model_pathways_to_cytoscope_test(models: &[GoCamModel])
   -> CytoscapeModelElements
{
    let mut model_map = HashMap::new();

    for model in models.iter() {
        model_map.insert(model.id().to_owned(), model);
    }

    let overlaps = GoCamModel::find_overlaps(&models);

    let mut edges = BTreeSet::new();
    let mut nodes = BTreeSet::new();

    let mut models_in_overlaps = HashSet::new();

    for overlap in &overlaps {
        let overlap_id = overlap.id();

        let overlap_node =
            CytoscapeModel {
                id: overlap_id.clone(),
                title: overlap.node_label.to_owned(),
            };

        nodes.insert(overlap_node);

        for (overlap_model_id, _, _) in overlap.models.iter() {

            models_in_overlaps.insert(overlap_model_id.to_owned());

            let edge = CytoscapeModelConnection {
                id: format!("{}-{}", overlap_model_id, overlap_id),
                label: overlap.node_label.to_owned(),
                source: overlap_model_id.to_owned(),
                target: overlap_id.to_owned(),
                enabler_id: "".to_owned(),
                has_direction: false,
            };
            edges.insert(edge);
        }
    }

    nodes.extend(models_in_overlaps.iter()
        .map(|model_id| {
            let model = model_map.get(model_id).unwrap();

            CytoscapeModel {
                id: model.id().to_owned(),
                title: model.title().to_owned(),
            }
        }));

    CytoscapeModelElements{
        nodes,
        edges,
    }
}

pub fn find_holes(model: &GoCamModel) -> Vec<GoCamNode> {
    let node_iter = model.node_iterator();
    node_iter.filter_map(|(_, node)| {
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
        let Some(individual_type) = individual.get_individual_type()
        else {
            continue;
        };

        if individual.individual_is_gene() {
            gene_map.insert(individual.id.clone(), individual_type);
        }
    }

    for fact in model.facts() {
        gene_map.remove(&fact.object);
        gene_map.remove(&fact.subject);
    }

    gene_map.iter().map(|(k, v)| { (k.to_owned(), v.id().to_owned(), v.label().to_owned()) }).collect()
}


pub type ChadoData = BTreeMap<ModelId, ChadoModelData>;

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct ChadoModelData {
    pub complex_terms: BTreeSet<String>,
    pub contributors: BTreeSet<String>,
    pub date: String,
    pub genes: BTreeSet<String>,
    pub modified_gene_pro_terms: BTreeSet<String>,
    pub process_terms: BTreeSet<String>,
    pub occurs_in_terms: BTreeSet<String>,
    pub located_in_terms: BTreeSet<String>,
    pub title: String,
    pub title_terms: BTreeSet<String>,
}

lazy_static! {
    static ref TERM_IN_TITLE_RE: Regex = Regex::new(r"\((GO:\d\d\d\d+)\)").unwrap();
}

fn chado_data_helper(model: &GoCamModel) -> ChadoModelData {
    let mut complex_terms = BTreeSet::new();
    let mut occurs_in_terms = BTreeSet::new();
    let mut located_in_terms = BTreeSet::new();
    let mut genes = BTreeSet::new();
    let mut modified_gene_pro_terms = BTreeSet::new();
    let mut process_terms = BTreeSet::new();

    let mut add_gene = |g: &str| genes.insert(g.replace("PomBase:", ""));

    for (_, node) in model.node_iterator() {
        let process_termids = node.part_of_process.iter()
            .map(|p| p.id.clone());
        process_terms.extend(process_termids);

        for occurs_in in &node.occurs_in {
            occurs_in_terms.insert(occurs_in.id().to_owned());
            if let GoCamComponent::ComplexComponent(it) = occurs_in {
                complex_terms.insert(it.id().to_owned());
            }
        }

        if let Some(ref located_in) = node.located_in {
            located_in_terms.insert(located_in.id().to_owned());
            if let GoCamComponent::ComplexComponent(it) = located_in {
                complex_terms.insert(it.id().to_owned());
            }
        }

        match &node.node_type {
            GoCamNodeType::Unknown => (),
            GoCamNodeType::Chemical => (),
            GoCamNodeType::UnknownMRNA => (),
            GoCamNodeType::MRNA(mrna) => {
                if let Some(no_suffix) = mrna.id.strip_suffix(|c: char| c.is_numeric()) {
                    if let Some(no_suffix) = no_suffix.strip_suffix('.') {
                        add_gene(no_suffix);
                    }
                }
            },
            GoCamNodeType::Gene(gene) => {
                add_gene(gene.id());
            },
            GoCamNodeType::ModifiedProtein(modified_protein_termid) => {
                modified_gene_pro_terms.insert(modified_protein_termid.id().to_owned());
            },
            GoCamNodeType::Activity(enabled_by) => match enabled_by {
                GoCamEnabledBy::Chemical(_) => (),
                GoCamEnabledBy::Gene(gene) => {
                    add_gene(gene.id());
                },
                GoCamEnabledBy::ModifiedProtein(modified_protein_termid) => {
                    modified_gene_pro_terms.insert(modified_protein_termid.id().to_owned());
                },
                GoCamEnabledBy::Complex(complex) => {
                    complex_terms.insert(complex.id().to_owned());
                    for gene in &complex.has_part_genes {
                        add_gene(&gene);
                    }
                },
            }
        }
    }

    let title = model.title().to_owned();
    let mut title_terms = BTreeSet::new();

    for (_, [title_termid]) in TERM_IN_TITLE_RE.captures_iter(&title).map(|c| c.extract()) {
        title_terms.insert(title_termid.to_owned());
        process_terms.insert(title_termid.to_owned());
    }

    let contributors = model.contributors().iter()
        .map(|c| c.replace("https://orcid.org/", "")).collect();

    ChadoModelData {
        title,
        title_terms,
        date: model.date().to_owned(),
        contributors,
        complex_terms,
        occurs_in_terms,
        located_in_terms,
        genes,
        modified_gene_pro_terms,
        process_terms,
    }
}

pub fn make_chado_data(models: &[GoCamModel]) -> ChadoData {
    let mut ret = BTreeMap::new();

    for model in models {
        let id = model.id().replace("gomodel:", "");
        ret.insert(id, chado_data_helper(model));
    }

    ret
}
