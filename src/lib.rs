use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

extern crate serde_json;
#[macro_use] extern crate serde_derive;

#[macro_use] extern crate lazy_static;

use itertools::Itertools;

use petgraph::{graph::NodeIndex, Undirected};
use petgraph::visit::Bfs;
use petgraph::visit::EdgeRef;

use pombase_gocam::{GoCamComponent, GoCamEnabledBy, GoCamModel, GoCamNode, GoCamNodeType, GoCamRawModel, ModelId};
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

pub fn get_connected_genes(model: &GoCamModel, min_connected_count: usize) -> HashSet<String> {
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
                ret.insert(gene.to_owned());
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

pub fn model_pathways_to_cytoscope(models: &[&GoCamModel]) -> CytoscapeElements {
    let mut model_map = HashMap::new();

    for model in models.iter() {
        model_map.insert(model.id().to_owned(), model);
    }

    let overlaps = GoCamModel::find_activity_overlaps(&models);

    let mut edges = vec![];

    let mut models_in_overlaps = HashSet::new();

    for overlap in &overlaps {
        let iter = overlap.model_ids.iter()
            .cartesian_product(overlap.model_ids.iter());

        for (first, second) in iter.into_iter() {
            if first >= second {
                continue;
            }
            models_in_overlaps.insert(first.to_owned());
            models_in_overlaps.insert(second.to_owned());
            let edge = CytoscapeEdge {
                data: CytoscapeEdgeData {
                    id: format!("{}-{}", first, second),
                    label: overlap.node_label.to_owned(),
                    source: first.to_owned(),
                    target: second.to_owned(),
                    weight: 0,
                }
            };
            edges.push(edge);
        }
    }

    let nodes: Vec<_> = models_in_overlaps.iter()
        .map(|model_id| {
            let model = model_map.get(model_id).unwrap();

            CytoscapeNode {
                data: CytoscapeNodeData {
                    id: model.id().to_owned(),
                    label: model.title().to_owned(),
                    db_id: "".to_owned(),
                    display_label: "".to_owned(),
                    type_string: "".to_owned(),
                    enabler_label: "".to_owned(),
                }
            }
        })
        .collect();

    CytoscapeElements {
        nodes,
        edges,
    }
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

    for node in model.node_iterator() {
        let process_termids = node.part_of_process.iter()
            .map(|i| i.id().to_owned());
        process_terms.extend(process_termids);

        if let Some(ref occurs_in) = node.occurs_in {
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

pub fn make_chado_data(models: &[&GoCamModel]) -> ChadoData {
    let mut ret = BTreeMap::new();

    for model in models {
        let id = model.id().replace("gomodel:", "");
        ret.insert(id, chado_data_helper(model));
    }

    ret
}
