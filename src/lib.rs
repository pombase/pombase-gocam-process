use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

extern crate serde_json;
#[macro_use] extern crate serde_derive;

#[macro_use] extern crate lazy_static;

use itertools::Itertools;

use petgraph::Direction;
use petgraph::{graph::NodeIndex, Undirected};
use petgraph::visit::Bfs;
use petgraph::visit::EdgeRef;

use pombase_gocam::overlaps::{find_activity_overlaps, GoCamNodeOverlap};
use pombase_gocam::{GoCamActivity, GoCamChemical, GoCamComplex, GoCamDirection, GoCamGeneIdentifier, GoCamGeneName, GoCamInput, GoCamModelTitle, GoCamOutput};
use pombase_gocam::{GoCamComponent, GoCamEnabledBy, GoCamModel,
                    GoCamNode, GoCamProcess,
                    GoCamNodeType, raw::GoCamRawModel, GoCamModelId};
use regex::Regex;

pub struct GoCamStats {
    pub id: GoCamModelId,
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
            if gocam_node.is_activity() {
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
                GoCamNodeType::Activity(GoCamActivity { ref enabler, inputs: _, outputs: _ }) => {
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
               .or_default()
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
    pub enabler_part_of_complex: Option<GoCamComplex>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub located_in: Option<GoCamComponent>,
    #[serde(skip_serializing_if="BTreeSet::is_empty")]
    pub occurs_in: BTreeSet<GoCamComponent>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub part_of_process: Option<GoCamProcess>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub happens_during: Option<GoCamProcess>,
    pub has_part_genes: BTreeSet<GoCamGeneIdentifier>,
    pub has_input: BTreeSet<GoCamInput>,
    pub has_output: BTreeSet<GoCamOutput>,
    // this node is in more than one model
    pub is_connecting_node: bool,
    #[serde(rename = "type")]
    pub type_string: String,
    #[serde(skip_serializing_if="Option::is_none")]
    pub parent: Option<GoCamModelId>,
    pub models: BTreeSet<(GoCamModelId, GoCamModelTitle)>,
}

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct CytoscapeEdgeData {
    pub id: CytoscapeId,
    pub label: String,
    pub source: CytoscapeId,
    pub target: CytoscapeId,
    pub ideal_edge_length: Option<u32>,
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
pub struct CytoscapeGeneInfo {
    pub name: Option<String>,
}

pub type CytoscapeGeneInfoMap = BTreeMap<GoCamGeneIdentifier, GoCamGeneName>;

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct CytoscapeElements {
    pub models: BTreeSet<(GoCamModelId, GoCamModelTitle)>,
    pub gene_info_map: CytoscapeGeneInfoMap,
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
    pub source: GoCamModelId,
    pub target: GoCamModelId,
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
                    ideal_edge_length: None,
                    weight: 0,
                }
            }
        }).collect();

    let nodes: Vec<_> = model.individuals()
        .filter_map(|individual| {
            if !seen_nodes.contains(&individual.id) {
                return None;
            }

            let individual_type = individual.types.first()?.to_owned();

            let (Some(ref label), Some(ref id)) = (individual_type.label, individual_type.id)
            else {
                return None;
            };

            let mut models = BTreeSet::new();
            models.insert((model.id().to_owned(), model.title()));

            Some(CytoscapeNode {
                data: CytoscapeNodeData {
                    id: individual.id.clone(),
                    type_string: "node".to_owned(),
                    display_label: label.to_owned(),
                    label: label.to_owned(),
                    node_id: id.clone(),
                    enabler_label: "".to_owned(),
                    enabler_id: "".to_owned(),
                    enabler_part_of_complex: None,
                    located_in: None,
                    occurs_in: BTreeSet::new(),
                    part_of_process: None,
                    happens_during: None,
                    has_part_genes: BTreeSet::new(),
                    has_input: BTreeSet::new(),
                    has_output: BTreeSet::new(),
                    is_connecting_node: false,
                    parent: None,
                    models,
                }
            })
        }).collect();

     CytoscapeElements {
        models: BTreeSet::new(),
        gene_info_map: BTreeMap::new(),
        nodes,
        edges,
    }
}

fn remove_suffix<'a>(s: &'a str, suffix: &str) -> &'a str {
    match s.strip_suffix(suffix) {
        Some(s) => s,
        None => s
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum GoCamCytoscapeStyle {
    IncludeParents,
    HideParents,
}

pub fn model_to_cytoscape_simple(model: &GoCamModel, overlaps: &Vec<GoCamNodeOverlap>,
                                 style: GoCamCytoscapeStyle) -> CytoscapeElements {
    let mut merged_models = BTreeSet::new();

    let mut model_map = BTreeMap::new();

    for overlap in overlaps {
        for individual_id in &overlap.overlapping_individual_ids {
            let overlap_models: BTreeSet<(GoCamModelId, GoCamModelTitle)>  = overlap.models.iter()
                .map(|(model_id, model_title, _)| {
                    (model_id.clone(), model_title.clone())
                })
                .collect();
           model_map.entry(individual_id.clone())
               .or_insert_with(BTreeSet::new)
               .extend(overlap_models.into_iter());
        }
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

            // allow longer edges between models
            let ideal_edge_length =
                if style == GoCamCytoscapeStyle::HideParents {
                    150
                } else if subject_node.original_model_id == object_node.original_model_id {
                    70
                } else {
                    350
                };

            CytoscapeEdge {
                data: CytoscapeEdgeData {
                    id: edge.fact_gocam_id.clone(),
                    label,
                    source,
                    target,
                    ideal_edge_length: Some(ideal_edge_length),
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
                } else if !enabler_label.is_empty() {
                    remove_suffix(remove_suffix(enabler_label, " Spom"), " Dmel").to_owned()
                } else {
                    "".to_owned()
                };
            let enabler_id = node.enabler_id().to_owned();
            let display_label =
                if !enabler_label.is_empty() {
                    enabler_label.clone()
                } else {
                    label.clone()
                };

            let node_models = model_map.get(node.source_ids.first().unwrap()).cloned()
                .unwrap_or_else(|| {
                    node.models.clone()
                });

            let is_connecting_node = node_models.len() > 1;

            let has_part_genes =
                if let GoCamNodeType::Activity(GoCamActivity { ref enabler, inputs: _, outputs: _ }) = node.node_type {
                    match enabler {
                        GoCamEnabledBy::Complex(complex) => {
                            complex.has_part_genes.iter()
                                .map(|id| {
                                    if let Some(new_id) = id.split(':').nth(1) {
                                        new_id.to_owned()
                                    } else {
                                        id.to_owned()
                                    }
                                })
                                .collect()
                        }
                        _ => BTreeSet::new(),
                    }
                } else {
                    BTreeSet::new()
                };

            merged_models.extend(node.models.iter().map(|v| v.to_owned()));

            let model_ids: BTreeSet<_> = node.models.iter()
                .map(|(model_id, _)| model_id.to_owned())
                .collect();

            let parent = if style == GoCamCytoscapeStyle::IncludeParents {
                if model_ids.len() == 1 {
                    model_ids.first().map(String::to_owned)
                } else {
                    node.original_model_id.as_ref()
                        .map(|original_model_id| original_model_id.to_owned())
                }
            } else {
                None
            };

            let mut has_input = BTreeSet::new();
            let mut has_output = BTreeSet::new();

            if let GoCamNodeType::Activity(GoCamActivity { enabler: _, ref inputs, ref outputs }) = node.node_type {
                has_input = inputs.clone();
                has_output = outputs.clone();
            }

            let enabler_part_of_complex =
                if let GoCamNodeType::Activity(GoCamActivity { ref enabler, inputs: _, outputs: _ }) = node.node_type &&
                   let GoCamEnabledBy::Gene(gene) = enabler &&
                   let Some(ref complex) = gene.part_of_complex
            {
                Some(complex.clone())
            } else {
                None
            };

            let located_in =
                if let GoCamNodeType::Chemical(ref chemical) = node.node_type {
                    chemical.located_in.clone()
                } else {
                    None
                };

            CytoscapeNode {
                data: CytoscapeNodeData {
                    id: node.individual_gocam_id.clone(),
                    type_string: node.type_string().to_owned(),
                    display_label,
                    label,
                    node_id: node.node_id.clone(),
                    enabler_label,
                    enabler_id,
                    enabler_part_of_complex,
                    located_in,
                    occurs_in: node.occurs_in.clone(),
                    part_of_process: node.part_of_process.clone(),
                    happens_during: node.happens_during.clone(),
                    has_part_genes,
                    has_input,
                    has_output,
                    is_connecting_node,
                    parent,
                    models: node_models.clone(),
                }
            }
        }).collect();

    let gene_info_map = model.genes_enabling_activities()
        .iter()
        .filter_map(|(id, name)| {
            name.as_ref().map(|name| (id.to_owned(), name.to_owned()))
        })
        .collect();


    CytoscapeElements {
        models: merged_models,
        gene_info_map,
        nodes,
        edges,
    }
}

pub fn model_connections_to_cytoscope(overlaps: &Vec<GoCamNodeOverlap>, model_ids_and_titles: &[(String, String)])
   -> CytoscapeModelElements
{
    let mut edges = BTreeSet::new();

    let mut models_in_overlaps = HashSet::new();

    for overlap in overlaps {
        let enabler_or_chemical_id =
            if let GoCamNodeType::Activity(GoCamActivity { ref enabler, inputs: _, outputs: _ }) = overlap.node_type {
                enabler.id().to_owned()
            } else if let GoCamNodeType::Chemical(ref chemical) = overlap.node_type {
                chemical.id().to_owned()
            } else {
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

/*
            eprintln!("{} {} {} {}",
                      first_model_id.to_string(), first_direction.to_string(),
                      second_model_id.to_string(), second_direction.to_string());
*/

            let (source, target) =
                if first_direction == second_direction {
                    continue;
                } else if (first_direction == GoCamDirection::IncomingConstitutivelyUpstream  ||
                    first_direction == GoCamDirection::Incoming ||
                    second_direction == GoCamDirection::Outgoing) &&
                   second_direction != GoCamDirection::IncomingConstitutivelyUpstream
                {
                    (first_model_id, second_model_id).to_owned()
                } else {
                    (second_model_id, first_model_id).to_owned()
                };

            let enabler_id = enabler_or_chemical_id.clone();

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

    let overlaps = find_activity_overlaps(models);

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

/// Return a Vec of nodes that have an unknown enabler
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

pub fn find_detached_chemicals(model: &GoCamModel) -> Vec<GoCamChemical> {
    let mut results = vec![];

    for (node_idx, node) in model.node_iterator() {
        if let GoCamNodeType::Chemical(ref chemical) = node.node_type {
            let in_iter = model.graph().edges_directed(node_idx, Direction::Incoming);
            let out_iter = model.graph().edges_directed(node_idx, Direction::Outgoing);
            let edge_iter = in_iter.chain(out_iter);

            if edge_iter.count() == 0 {
                results.push(chemical.clone());
            }
        }
    }

    results
}


pub type ChadoData = BTreeMap<GoCamModelId, ChadoModelData>;

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct ChadoModelData {
    pub complex_terms: BTreeSet<String>,
    pub contributors: BTreeSet<String>,
    pub date: String,
    pub genes: BTreeSet<String>,
    pub target_genes: BTreeSet<String>,
    pub modified_gene_pro_terms: BTreeSet<String>,
    pub modified_target_gene_pro_terms: BTreeSet<String>,
    pub process_terms: BTreeSet<String>,
    pub occurs_in_terms: BTreeSet<String>,
    pub located_in_terms: BTreeSet<String>,
    pub title: String,
    pub title_terms: BTreeSet<String>,
    pub pathway_holes: Vec<GoCamNode>,
}

lazy_static! {
    static ref TERM_IN_TITLE_RE: Regex = Regex::new(r"\(\s*(GO:\d\d\d\d+)\s*\)").unwrap();
}

pub fn chado_data_helper(model: &GoCamModel) -> ChadoModelData {
    let mut complex_terms = BTreeSet::new();
    let mut occurs_in_terms = BTreeSet::new();
    let mut located_in_terms = BTreeSet::new();
    let mut target_genes = BTreeSet::new();
    let mut modified_gene_pro_terms = BTreeSet::new();
    let mut modified_target_gene_pro_terms = BTreeSet::new();
    let mut process_terms = BTreeSet::new();

    let mut add_target = |g: &str| target_genes.insert(g.replace("PomBase:", ""));

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

        if let GoCamNodeType::Chemical(ref chemical) = node.node_type {
            if let Some(ref located_in) = chemical.located_in {
                located_in_terms.insert(located_in.id().to_owned());
                if let GoCamComponent::ComplexComponent(it) = located_in {
                    complex_terms.insert(it.id().to_owned());
                }
            }
        }

        match &node.node_type {
            GoCamNodeType::Unknown => (),
            GoCamNodeType::Chemical(_) => (),
            GoCamNodeType::UnknownMRNA => (),
            GoCamNodeType::MRNA(mrna) => {
                if let Some(no_suffix) = mrna.id.strip_suffix(|c: char| c.is_numeric()) {
                    if let Some(no_suffix) = no_suffix.strip_suffix('.') {
                        add_target(no_suffix);
                    }
                }
            },
            GoCamNodeType::Gene(gene) => {
                add_target(gene.id());
            },
            GoCamNodeType::ModifiedProtein(modified_protein_termid) => {
                modified_target_gene_pro_terms.insert(modified_protein_termid.id().to_owned());
            },
            GoCamNodeType::Activity(GoCamActivity { enabler, inputs: _, outputs: _ }) => match enabler {
                GoCamEnabledBy::Chemical(_) => (),
                GoCamEnabledBy::Gene(_) => (),
                GoCamEnabledBy::ModifiedProtein(modified_protein_termid) => {
                    modified_gene_pro_terms.insert(modified_protein_termid.id().to_owned());
                },
                GoCamEnabledBy::Complex(complex) => {
                    complex_terms.insert(complex.id().to_owned());
                },
            },
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

    let pathway_holes = find_holes(model);

    let genes = model.genes_enabling_activities().keys()
        .map(|g| g.replace("PomBase:", "")).collect();

    ChadoModelData {
        title,
        title_terms,
        date: model.date().to_owned(),
        contributors,
        complex_terms,
        occurs_in_terms,
        located_in_terms,
        genes,
        target_genes,
        modified_gene_pro_terms,
        modified_target_gene_pro_terms,
        process_terms,
        pathway_holes,
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

#[cfg(test)]
mod tests {
    use std::{collections::HashMap, fs::File};

    use pombase_gocam::{parse_gocam_model, GoCamActivity, GoCamEnabledBy, GoCamNodeType};

    use crate::{chado_data_helper, find_holes};

    #[test]
    fn find_holes_test() {
        let mut source = File::open("tests/data/gomodel_67c10cc400002026.json").unwrap();
        let model = parse_gocam_model(&mut source).unwrap();
        assert_eq!(model.id(), "gomodel:67c10cc400002026");

        let holes = find_holes(&model);

        assert_eq!(holes.len(), 1);
        let first_hole = holes.first().unwrap();
        let GoCamNodeType::Activity(GoCamActivity { enabler, inputs, outputs }) = &first_hole.node_type
        else {
            panic!();
        };

        if let GoCamEnabledBy::Chemical(chemical) = enabler {
            assert_eq!(chemical.id(), "CHEBI:36080");
        } else {
            panic!();
        }

        assert_eq!(inputs.len(), 1);
        assert_eq!(outputs.len(), 1);

        let input = inputs.first().unwrap();
        assert_eq!(input.id(), "CHEBI:149473");
        assert_eq!(input.located_in.clone().unwrap().id(),
                   "GO:0005829");
        let output = outputs.first().unwrap();
        assert_eq!(output.id(), "CHEBI:149473");
        assert_eq!(output.located_in.clone().unwrap().label(),
                   "mitochondrion");
    }

    #[test]
    fn test_chado_data_helper() {
        let mut source = File::open("tests/data/gomodel_66187e4700001744.json").unwrap();
        let mut model = parse_gocam_model(&mut source).unwrap();

        let mut pro_term_to_gene_map = HashMap::new();

        pro_term_to_gene_map.insert("PR:000059631".to_owned(),
                                    "PomBase:SPCC1020.02".to_owned());
        pro_term_to_gene_map.insert("PR:000027566".to_owned(),
                                    "PomBase:SPCC622.08c".to_owned());
        pro_term_to_gene_map.insert("PR:000027557".to_owned(),
                                    "PomBase:SPAC19G12.06c".to_owned());
        pro_term_to_gene_map.insert("PR:000050512".to_owned(),
                                    "PomBase:SPBC29A10.14".to_owned());

        model.add_pro_term_to_gene_map(&pro_term_to_gene_map);

        assert_eq!(model.id(), "gomodel:66187e4700001744");

        let expected_genes_in_model = vec![
            "SPAC15E1.07c".to_owned(),
            "SPAC19G12.06c".to_owned(),
            "SPAC23C11.16".to_owned(),
            "SPAC664.01c".to_owned(),
            "SPBC106.01".to_owned(),
            "SPBC16H5.07c".to_owned(),
            "SPBC29A10.14".to_owned(),
            "SPBP35G2.03c".to_owned(),
            "SPCC1020.02".to_owned(),
            "SPCC1322.12c".to_owned(),
            "SPCC622.08c".to_owned(),
        ];

        let chado_data = chado_data_helper(&model);

        let mut chado_data_genes: Vec<_> = chado_data.genes.iter().cloned().collect();

        chado_data_genes.sort();

        assert_eq!(expected_genes_in_model, chado_data_genes);
    }

}
