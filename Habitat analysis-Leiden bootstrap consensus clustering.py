import warnings
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from numpy.random import default_rng
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

warnings.filterwarnings("ignore")

try:
    import igraph as ig
    import leidenalg
except ImportError as exc:  # pragma: no cover - runtime dependency check
    raise ImportError(
        "Leiden clustering requires the `igraph` and `leidenalg` packages.\n"
        "Install them before running this script, e.g.:\n"
        "    pip install python-igraph leidenalg"
    ) from exc


METADATA_COLUMNS = [
    "Center",
    "patient",
    "superpixel",
    "superpixel_id",
    "superpixel_id_original",
    "Sex",
    "Age",
    "BMI",
    "Smoking",
    "Disease",
    "PET_SUL_mean",
    "PET_Entropy_mean",
    "CT_HU_mean",
    "CT_Entropy_mean",
]

DISEASE_LABELS = {
    1: "Benign",
    2: "ADC",
    3: "SqCC",
    4: "SCLC",
}

SMOKING_LABELS = {
    0: "Non-smoker",
    1: "Smoker",
}

DEFAULT_FEATURE_COLUMNS = [
    "PET_SUL_mean",
    "PET_Entropy_mean",
    "CT_HU_mean",
    "CT_Entropy_mean",
]


def load_superpixel_features(feature_path: Path) -> pd.DataFrame:
    """
    Load the GLM-adjusted superpixel-level feature matrix.

    Parameters
    ----------
    feature_path : Path
        Excel file that stores superpixel metadata and features.

    Returns
    -------
    pd.DataFrame
        Feature table where each row corresponds to a superpixel.
    """

    print("\n" + "=" * 80)
    print("加载超像素特征数据")
    print("=" * 80)

    df = pd.read_excel(feature_path)
    print(f"总超像素数: {len(df)}")

    # Normalise metadata columns to numeric types where possible
    for col in ["Disease", "Smoking"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")

    metadata_preview = [col for col in METADATA_COLUMNS if col in df.columns][:10]
    print(f"元数据列示例: {metadata_preview}")

    n_patients = df["patient"].nunique()
    print(f"患者数量: {n_patients}")

    disease_counts = df.dropna(subset=["Disease"])["Disease"].value_counts()
    print("疾病类型分布 (按超像素计数):")
    for code, count in disease_counts.sort_index().items():
        name = DISEASE_LABELS.get(int(code), f"Disease {int(code)}")
        print(f"  {name} ({int(code)}): {count}")

    return df


def select_balanced_superpixels(
    df: pd.DataFrame,
    include_all_diseases: Iterable[int] = (1, 4),
    sample_diseases: Iterable[int] = (2, 3),
    n_per_smoking: int = 6,
    random_state: int = 42,
) -> Tuple[pd.DataFrame, Dict[str, Dict[str, int]]]:
    """
    Select a balanced set of patients/superpixels according to the specification.

    Rules
    -----
    * Include all superpixels from patients whose Disease is in `include_all_diseases`.
    * For each disease in `sample_diseases`, randomly choose `n_per_smoking` patients
      for Smoking=0 and for Smoking=1 (if available), and include all of their superpixels.

    Parameters
    ----------
    df : pd.DataFrame
        Superpixel feature table.
    include_all_diseases : Iterable[int]
        Diseases for which all patients are retained.
    sample_diseases : Iterable[int]
        Diseases for which patients are subsampled per smoking group.
    n_per_smoking : int
        Number of patients to sample for each smoking state per disease.
    random_state : int
        Random seed for reproducibility.

    Returns
    -------
    Tuple[pd.DataFrame, Dict[str, Dict[str, int]]]
        Filtered DataFrame containing the selected superpixels, and a summary
        dictionary describing how many patients were included per group.
    """

    rng = default_rng(random_state)
    selected_patients = set()
    selection_summary: Dict[str, Dict[str, int]] = {}

    df_meta = df[["patient", "Disease", "Smoking"]].drop_duplicates()

    # Include all patients from specified diseases
    include_all_diseases = tuple(include_all_diseases)
    include_mask = df_meta["Disease"].isin(include_all_diseases)
    include_patients = df_meta.loc[include_mask, "patient"].unique().tolist()
    selected_patients.update(include_patients)
    for disease in include_all_diseases:
        key = f"Disease_{disease}"
        selection_summary.setdefault(key, {})
        n_patients = df_meta.loc[df_meta["Disease"] == disease, "patient"].nunique()
        selection_summary[key]["all"] = int(n_patients)

    # Sample per smoking status for the remaining diseases
    sample_diseases = tuple(sample_diseases)
    for disease in sample_diseases:
        disease_key = f"Disease_{disease}"
        selection_summary.setdefault(disease_key, {})
        for smoking in [0, 1]:
            mask = (df_meta["Disease"] == disease) & (df_meta["Smoking"] == smoking)
            candidates = df_meta.loc[mask, "patient"].unique().tolist()
            if len(candidates) == 0:
                selection_summary[disease_key][f"Smoking_{smoking}"] = 0
                print(
                    f"WARNING: 无可用患者 (Disease={disease}, Smoking={smoking})，"
                    "该组合将在聚类中缺失。"
                )
                continue

            n_select = min(n_per_smoking, len(candidates))
            chosen = rng.choice(candidates, size=n_select, replace=False)
            selected_patients.update(chosen)
            selection_summary[disease_key][f"Smoking_{smoking}"] = int(n_select)

    df_selected = df[df["patient"].isin(selected_patients)].reset_index(drop=True)
    print("\n" + "=" * 80)
    print("选择用于聚类的超像素")
    print("=" * 80)
    print(f"选择患者数: {df_selected['patient'].nunique()}")
    print(f"选择超像素数: {len(df_selected)}")
    print("选择摘要:")
    for disease_key, stats in selection_summary.items():
        print(f"  {disease_key}: {stats}")

    return df_selected, selection_summary


def prepare_superpixel_matrix(
    df_selected: pd.DataFrame,
) -> Tuple[pd.DataFrame, np.ndarray, List[str]]:
    """
    Prepare the superpixel feature matrix for clustering.

    Parameters
    ----------
    df_selected : pd.DataFrame
        Filtered DataFrame containing only the selected superpixels.

    Returns
    -------
    Tuple[pd.DataFrame, np.ndarray, List[str]]
        - Metadata DataFrame (one row per superpixel) with a unique superpixel_id column.
        - Scaled feature matrix (numpy array).
        - List of feature column names used for clustering.
    """

    df_meta = df_selected.copy()

    if "superpixel_id" in df_meta.columns:
        df_meta = df_meta.rename(columns={"superpixel_id": "superpixel_id_original"})
    else:
        df_meta["superpixel_id_original"] = (
            df_meta["superpixel"].astype(str)
            if "superpixel" in df_meta.columns
            else np.nan
        )

    # Create deterministic, patient-specific unique IDs to avoid collisions.
    seq = df_meta.groupby("patient").cumcount().add(1)
    df_meta["superpixel_id"] = (
        df_meta["patient"].astype(str)
        + "_SP"
        + seq.astype(str)
    )

    df_meta["superpixel_id"] = df_meta["superpixel_id"].astype(str)

    feature_cols = [
        col for col in DEFAULT_FEATURE_COLUMNS if col in df_meta.columns
    ]

    if not feature_cols:
        feature_cols = [
            col
            for col in df_meta.columns
            if col not in METADATA_COLUMNS + ["superpixel_id"]
        ]

    if len(feature_cols) == 0:
        raise ValueError(
            "未找到用于聚类的特征列。请确认数据包含 PET/CT 特征或移除它们在元数据配置中的声明。"
        )

    scaler = StandardScaler()
    features_scaled = scaler.fit_transform(df_meta[feature_cols].values)

    print("\n" + "=" * 80)
    print("准备超像素特征矩阵")
    print("=" * 80)
    print(f"特征数量: {len(feature_cols)}")
    print(f"矩阵形状: {features_scaled.shape}")

    return df_meta, features_scaled, feature_cols


def determine_consensus_threshold(
    consensus_matrix: np.ndarray,
    *,
    base_threshold: float = 0.2,
    quantile: float = 0.7,
) -> float:
    """
    Determine an adaptive edge保留阈值 based on the consensus matrix distribution.
    """

    if consensus_matrix.size == 0:
        return base_threshold

    upper_idx = np.triu_indices_from(consensus_matrix, k=1)
    weights = consensus_matrix[upper_idx]
    weights = weights[np.isfinite(weights)]
    weights = weights[weights > 0]

    if len(weights) == 0:
        return base_threshold

    quantile_value = float(np.quantile(weights, quantile))
    threshold = max(base_threshold, quantile_value)
    return threshold


def build_similarity_graph(
    similarity_matrix: np.ndarray,
    node_ids: List[str],
    k_neighbors: int = 15,
    min_weight: float = 0.0,
    verbose: bool = True,
    k_cross: Optional[int] = None,
    use_jaccard: bool = True,
) -> nx.Graph:
    """
    Construct a mutual k-nearest-neighbour graph. Optionally uses Jaccard overlap
    for edge weights; otherwise retains the original similarity values.

    Parameters
    ----------
    similarity_matrix : np.ndarray
        Pairwise similarity / affinity matrix.
    node_ids : list[str]
        Identifiers that map each row/column to a node label.
    k_neighbors : int
        Number of nearest neighbours considered per node before applying the
        mutual filter.
    min_weight : float
        Drop edges whose Jaccard weight (shared-top-k overlap) is below this
        threshold.
    verbose : bool
        Whether to print diagnostic information.
    k_cross : int or None
        Size of the neighbour set used to evaluate the Jaccard overlap. If None,
        Jaccard weighting is skipped.
    use_jaccard : bool
        Whether to compute weights via Jaccard overlap. If False, the original
        similarity values are used.
    """

    n_nodes = len(node_ids)
    if similarity_matrix.shape != (n_nodes, n_nodes):
        raise ValueError("similarity_matrix 与 node_ids 长度不匹配。")

    G = nx.Graph()
    for node in node_ids:
        G.add_node(node)

    if n_nodes <= 1:
        return G

    k = max(1, min(k_neighbors, n_nodes - 1))
    if k_cross is not None:
        k_cross = max(1, min(k_cross, k))

    top_indices_per_node: List[List[int]] = []
    cross_sets: List[set] = []

    for i in range(n_nodes):
        sims = similarity_matrix[i].copy()
        sims[i] = -np.inf  # exclude self
        top_indices = np.argsort(sims)[-k:]
        top_indices_per_node.append(top_indices.tolist())

        if k_cross is not None:
            cross_neighbors = top_indices[-k_cross:].tolist()
            cross_sets.append(set(cross_neighbors))

    for i in range(n_nodes):
        for j in top_indices_per_node[i]:
            if i == j:
                continue
            # mutual kNN check
            if i not in top_indices_per_node[j]:
                continue

            if use_jaccard and k_cross is not None:
                neigh_i = cross_sets[i]
                neigh_j = cross_sets[j]
                union = neigh_i | neigh_j
                if not union:
                    continue
                weight = len(neigh_i & neigh_j) / len(union)
            else:
                weight = float(similarity_matrix[i, j])
            if weight < min_weight:
                continue

            u, v = node_ids[i], node_ids[j]
            if G.has_edge(u, v):
                G[u][v]["weight"] = max(G[u][v]["weight"], weight)
            else:
                G.add_edge(u, v, weight=weight)

    if verbose:
        print("\n" + "=" * 80)
        if use_jaccard and k_cross is not None:
            print("构建相似性图 (mutual kNN + Jaccard 权重)")
            print(f"k = {k}, k_cross = {k_cross}, min_weight = {min_weight:.2f}")
        else:
            print("构建相似性图 (mutual kNN + 原始相似度权重)")
            print(f"k = {k}, min_weight = {min_weight:.2f}")
        print("=" * 80)
        print(f"节点数: {G.number_of_nodes()}")
        print(f"边数: {G.number_of_edges()}")
        print(f"k = {k}, k_cross = {k_cross}, min_weight = {min_weight:.2f}")

    return G


def perform_leiden_clustering(
    G: nx.Graph,
    resolution: float = 0.8,
    random_state: int = 42,
    verbose: bool = True,
) -> Tuple[Dict[str, int], float]:
    """
    Run Leiden community detection on a NetworkX graph.
    """

    if G.number_of_nodes() == 0:
        return {}, 0.0

    nodes = list(G.nodes())
    node_to_idx = {node: idx for idx, node in enumerate(nodes)}

    ig_graph = ig.Graph()
    ig_graph.add_vertices(len(nodes))
    ig_graph.vs["name"] = nodes

    if G.number_of_edges() > 0:
        edges = [(node_to_idx[u], node_to_idx[v]) for u, v in G.edges()]
        weights = [float(G[u][v].get("weight", 1.0)) for u, v in G.edges()]
        ig_graph.add_edges(edges)
        ig_graph.es["weight"] = weights

    partition = leidenalg.find_partition(
        ig_graph,
        leidenalg.RBConfigurationVertexPartition,
        weights="weight" if G.number_of_edges() > 0 else None,
        resolution_parameter=resolution,
        seed=random_state,
    )

    membership = partition.membership
    modularity = getattr(partition, "modularity", 0.0)
    communities = {nodes[idx]: int(comm) for idx, comm in enumerate(membership)}

    if verbose:
        print("\n" + "=" * 80)
        print("执行 Leiden 社区检测")
        print("=" * 80)
        print(f"社区数量: {len(set(membership))}")
        print(f"模块度: {modularity:.4f}")

    return communities, modularity


def normalize_partition_labels(partition: Dict[str, int]) -> Dict[str, int]:
    """
    Relabel partition identifiers to a contiguous range starting at 0.
    """

    label_map = {
        label: idx for idx, label in enumerate(sorted(set(partition.values())))
    }
    return {node: label_map[label] for node, label in partition.items()}


def compute_partition_modularity(
    G: nx.Graph,
    partition: Dict[str, int],
) -> float:
    """
    Compute modularity for a given partition on the provided graph.
    """

    if G.number_of_nodes() == 0:
        return 0.0

    nodes = list(G.nodes())
    node_to_idx = {node: idx for idx, node in enumerate(nodes)}

    ig_graph = ig.Graph()
    ig_graph.add_vertices(len(nodes))
    ig_graph.vs["name"] = nodes

    weights: Optional[List[float]] = None
    if G.number_of_edges() > 0:
        edges = [(node_to_idx[u], node_to_idx[v]) for u, v in G.edges()]
        weights = [float(G[u][v].get("weight", 1.0)) for u, v in G.edges()]
        ig_graph.add_edges(edges)
        ig_graph.es["weight"] = weights

    membership = [partition[node] for node in nodes]
    if weights is not None:
        modularity = ig_graph.modularity(membership, weights=weights)
    else:
        modularity = ig_graph.modularity(membership)

    return float(modularity)


def merge_partition_to_target(
    partition: Dict[str, int],
    superpixel_ids: List[str],
    consensus_matrix: np.ndarray,
    target_k: int,
    method: str = "average",
) -> Tuple[Dict[str, int], bool]:
    """
    Merge existing communities using层次聚类 until reaching target_k clusters.
    """

    if target_k <= 0:
        raise ValueError("target_k must be positive.")

    unique_labels = sorted(set(partition.values()))
    if len(unique_labels) <= target_k:
        return partition, False

    id_to_idx = {sp_id: idx for idx, sp_id in enumerate(superpixel_ids)}
    clusters: List[int] = []
    member_indices: List[List[int]] = []
    for label in unique_labels:
        indices = [id_to_idx[sp_id] for sp_id in superpixel_ids if partition[sp_id] == label]  # noqa: E501
        if not indices:
            continue
        clusters.append(label)
        member_indices.append(indices)

    n_clusters = len(clusters)
    if n_clusters <= target_k:
        return partition, False

    cluster_sim = np.zeros((n_clusters, n_clusters), dtype=float)
    for i in range(n_clusters):
        idx_i = member_indices[i]
        cluster_sim[i, i] = float(consensus_matrix[np.ix_(idx_i, idx_i)].mean())
        for j in range(i + 1, n_clusters):
            idx_j = member_indices[j]
            sim = float(consensus_matrix[np.ix_(idx_i, idx_j)].mean())
            cluster_sim[i, j] = cluster_sim[j, i] = sim

    dist_matrix = 1.0 - cluster_sim
    dist_matrix = np.clip(dist_matrix, 0.0, 1.0)
    np.fill_diagonal(dist_matrix, 0.0)

    condensed = squareform(dist_matrix, checks=False)
    if condensed.size == 0:
        return partition, False

    linkage_matrix = linkage(condensed, method=method)
    merged_labels = fcluster(linkage_matrix, t=target_k, criterion="maxclust")

    label_map = {
        clusters[i]: int(merged_labels[i] - 1) for i in range(len(clusters))
    }
    merged_partition = {
        sp_id: label_map.get(partition[sp_id], partition[sp_id])
        for sp_id in superpixel_ids
    }

    merged_partition = normalize_partition_labels(merged_partition)
    return merged_partition, True


def search_leiden_resolutions(
    G: nx.Graph,
    resolutions: Iterable[float],
    target_range: Tuple[int, int],
    random_state: int = 42,
) -> Tuple[Optional[Dict[str, int]], Optional[float], Optional[float], Optional[int], pd.DataFrame]:
    """
    Evaluate Leiden clustering across a set of resolution values and pick the
    configuration that best matches the desired community count range.
    """

    if G.number_of_nodes() == 0:
        return None, None, None, None, pd.DataFrame()

    records: List[Dict[str, float]] = []
    best_partition: Optional[Dict[str, int]] = None
    best_modularity: Optional[float] = None
    best_resolution: Optional[float] = None
    best_n_communities: Optional[int] = None
    best_priority: Optional[Tuple[int, float]] = None

    resolutions = list(resolutions)
    lower, upper = target_range

    for idx, gamma in enumerate(resolutions):
        partition, modularity = perform_leiden_clustering(
            G,
            resolution=float(gamma),
            random_state=random_state + idx,
            verbose=False,
        )
        n_communities = len(set(partition.values()))
        diff = 0
        if not (lower <= n_communities <= upper):
            diff = min(abs(n_communities - lower), abs(n_communities - upper))

        priority = (diff, -modularity)
        if best_priority is None or priority < best_priority:
            best_priority = priority
            best_partition = partition
            best_modularity = modularity
            best_resolution = float(gamma)
            best_n_communities = n_communities

        records.append(
            {
                "resolution": float(gamma),
                "n_communities": n_communities,
                "modularity": modularity,
                "diff_from_target": diff,
            }
        )

    if best_partition is not None:
        best_partition = normalize_partition_labels(best_partition)
        best_n_communities = len(set(best_partition.values()))

    results_df = pd.DataFrame(records)
    return (
        best_partition,
        best_modularity,
        best_resolution,
        best_n_communities,
        results_df,
    )


def bootstrap_consensus_leiden(
    features: np.ndarray,
    node_ids: List[str],
    *,
    n_bootstrap: int = 300,
    sample_ratio: float = 0.8,
    k_neighbors: int = 15,
    resolution: float = 1.0,
    k_range: Optional[Tuple[int, int]] = None,
    resolution_range: Optional[Tuple[float, float]] = None,
    random_state: int = 42,
) -> Dict[str, np.ndarray]:
    """
    Build a consensus matrix through bootstrap resampling and Leiden clustering.

    Each bootstrap draws an 80% subset *without* replacement, perturbs the graph
    construction with a random k-neighbour value (within `k_range`) and a random
    Leiden resolution (`resolution_range`), then accumulates how often node pairs
    co-occur in the same community. This setup evaluates the stability of the
    clustering under both sampling and parameter jitter.
    """

    if not 0 < sample_ratio <= 1:
        raise ValueError("sample_ratio 必须位于 (0, 1] 区间内。")

    if k_range is not None:
        k_min, k_max = k_range
        if k_min <= 0 or k_max < k_min:
            raise ValueError("k_range 必须满足 0 < k_min <= k_max。")

    if resolution_range is not None:
        gamma_min, gamma_max = resolution_range
        if gamma_min <= 0 or gamma_max < gamma_min:
            raise ValueError("resolution_range 必须满足 0 < γ_min <= γ_max。")

    rng = default_rng(random_state)
    n_nodes = len(node_ids)
    node_to_pos = {node: idx for idx, node in enumerate(node_ids)}

    cooccurrence = np.zeros((n_nodes, n_nodes), dtype=float)
    sample_counts = np.zeros((n_nodes, n_nodes), dtype=float)
    modularities: List[float] = []

    min_subset = max(2, int(np.ceil(sample_ratio * n_nodes)))
    progress_interval = max(1, n_bootstrap // 10)

    print("\n" + "=" * 80)
    print("启动 bootstrap 共识聚类")
    print("=" * 80)
    print(f"Bootstrap 次数: {n_bootstrap}")
    print(f"每次抽样比例: {sample_ratio:.2f}")
    print(f"最小抽样数量: {min_subset}")

    k_history: List[int] = []
    gamma_history: List[float] = []

    for b in range(n_bootstrap):
        subset_indices = np.sort(
            rng.choice(n_nodes, size=min_subset, replace=False)
        )
        subset_ids = [node_ids[idx] for idx in subset_indices]
        subset_features = features[subset_indices]

        if k_range is not None:
            k_boot = int(rng.integers(k_range[0], k_range[1] + 1))
        else:
            k_boot = k_neighbors
        k_history.append(k_boot)

        if resolution_range is not None:
            resolution_boot = float(rng.uniform(resolution_range[0], resolution_range[1]))
        else:
            resolution_boot = resolution
        gamma_history.append(resolution_boot)

        sim_matrix = cosine_similarity(subset_features)
        G_subset = build_similarity_graph(
            sim_matrix,
            subset_ids,
            k_neighbors=k_boot,
            min_weight=0.2,
            k_cross=None,
            use_jaccard=False,
            verbose=False,
        )
        partition, modularity = perform_leiden_clustering(
            G_subset,
            resolution=resolution_boot,
            random_state=int(rng.integers(1_000_000_000)),
            verbose=False,
        )
        modularities.append(modularity)

        # Update sampling counts
        for i_pos, i_global in enumerate(subset_indices):
            sample_counts[i_global, i_global] += 1
            for j_pos in range(i_pos + 1, len(subset_indices)):
                j_global = subset_indices[j_pos]
                sample_counts[i_global, j_global] += 1
                sample_counts[j_global, i_global] += 1

        # Update co-occurrence counts
        communities: Dict[int, List[str]] = {}
        for node, comm in partition.items():
            communities.setdefault(comm, []).append(node)

        for members in communities.values():
            member_positions = [node_to_pos[m] for m in members]
            for idx in member_positions:
                cooccurrence[idx, idx] += 1
            for i_pos, i_global in enumerate(member_positions):
                for j_pos in range(i_pos + 1, len(member_positions)):
                    j_global = member_positions[j_pos]
                    cooccurrence[i_global, j_global] += 1
                    cooccurrence[j_global, i_global] += 1

        if (b + 1) % progress_interval == 0 or b == n_bootstrap - 1:
            print(
                f"  已完成 {b + 1}/{n_bootstrap} 次重采样，"
                f"k={k_boot}, γ={resolution_boot:.2f}, 模块度 {modularity:.4f}"
            )

    with np.errstate(divide="ignore", invalid="ignore"):
        consensus = np.divide(
            cooccurrence,
            sample_counts,
            out=np.zeros_like(cooccurrence),
            where=sample_counts > 0,
        )

    np.fill_diagonal(consensus, 1.0)

    return {
        "consensus": consensus,
        "cooccurrence": cooccurrence,
        "sample_counts": sample_counts,
        "modularities": np.array(modularities),
        "k_history": np.array(k_history, dtype=int),
        "resolution_history": np.array(gamma_history, dtype=float),
    }


def save_consensus_results(
    metadata: pd.DataFrame,
    partition: Dict[str, int],
    consensus_matrix: np.ndarray,
    sample_counts: np.ndarray,
    output_dir: Path,
) -> Dict[str, Path]:
    """
    Persist consensus clustering outputs and summary tables.
    """

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    metadata = metadata.copy()
    metadata["Community"] = metadata["superpixel_id"].map(partition)
    metadata["Community"] = metadata["Community"].astype(int) + 1

    assignments_path = output_dir / "leiden_consensus_superpixels.xlsx"
    metadata.to_excel(assignments_path, index=False)

    # Summary by community
    cluster_summary = metadata.groupby("Community").agg(
        N_Superpixels=("superpixel_id", "count"),
        N_Patients=("patient", "nunique"),
    )

    for disease_code in sorted(metadata["Disease"].dropna().unique()):
        disease_code = int(disease_code)
        column = f"Disease_{disease_code}"
        counts = (
            metadata[metadata["Disease"] == disease_code]
            .groupby("Community")["superpixel_id"]
            .count()
        )
        cluster_summary[column] = counts

    for smoking_code in sorted(metadata["Smoking"].dropna().unique()):
        smoking_code = int(smoking_code)
        column = f"Smoking_{smoking_code}"
        counts = (
            metadata[metadata["Smoking"] == smoking_code]
            .groupby("Community")["superpixel_id"]
            .count()
        )
        cluster_summary[column] = counts

    # Cross counts: disease x smoking
    cross_counts = (
        metadata.groupby(["Community", "Disease", "Smoking"])["superpixel_id"]
        .count()
        .rename("count")
        .reset_index()
    )
    for _, row in cross_counts.iterrows():
        disease_code = int(row["Disease"])
        smoking_code = int(row["Smoking"])
        col_name = f"D{disease_code}_S{smoking_code}"
        cluster_summary.loc[row["Community"], col_name] = row["count"]

    cluster_summary = cluster_summary.fillna(0).astype(int)
    summary_path = output_dir / "leiden_consensus_summary.xlsx"
    cluster_summary.to_excel(summary_path)

    # Consensus matrices
    superpixel_ids = metadata["superpixel_id"].tolist()
    consensus_df = pd.DataFrame(
        consensus_matrix,
        index=superpixel_ids,
        columns=superpixel_ids,
    )
    support_df = pd.DataFrame(
        sample_counts,
        index=superpixel_ids,
        columns=superpixel_ids,
    )

    consensus_path = output_dir / "consensus_matrix.csv"
    support_path = output_dir / "consensus_support_counts.csv"
    consensus_df.to_csv(consensus_path)
    support_df.to_csv(support_path)

    print("\n" + "=" * 80)
    print("结果已保存")
    print("=" * 80)
    print(f"  聚类分配: {assignments_path}")
    print(f"  社区汇总: {summary_path}")
    print(f"  共识矩阵: {consensus_path}")
    print(f"  抽样次数矩阵: {support_path}")

    return {
        "assignments": assignments_path,
        "summary": summary_path,
        "consensus_matrix": consensus_path,
        "support_matrix": support_path,
    }


def visualize_consensus_heatmap(
    consensus_matrix: np.ndarray,
    metadata: pd.DataFrame,
    partition: Dict[str, int],
    output_dir: Path,
) -> None:
    """
    Visualise the consensus matrix ordered by community assignment.
    """

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    community_series = metadata["superpixel_id"].map(partition)
    order = np.argsort(community_series.to_numpy())
    ordered_matrix = consensus_matrix[np.ix_(order, order)]

    plt.figure(figsize=(12, 10))
    sns.heatmap(
        ordered_matrix,
        cmap="viridis",
        vmin=0.0,
        vmax=1.0,
        cbar_kws={"label": "Consensus probability"},
    )
    plt.title("Leiden Bootstrap Consensus Heatmap", fontsize=16, fontweight="bold")
    plt.xlabel("Superpixel")
    plt.ylabel("Superpixel")
    plt.xticks([])
    plt.yticks([])
    plt.tight_layout()

    heatmap_path = output_dir / "consensus_heatmap.png"
    plt.savefig(heatmap_path, dpi=300)
    plt.close()

    print(f"  共识热图已保存: {heatmap_path}")


def visualize_network(
    G: nx.Graph,
    partition: Dict[str, int],
    output_dir: Path,
    max_nodes: int = 300,
) -> None:
    """
    Visualise the similarity network. For dense graphs only the top-`max_nodes`
    nodes (by degree) are displayed to keep the plot readable.
    """

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if G.number_of_nodes() == 0 or G.number_of_edges() == 0:
        print("  图结构为空，跳过网络可视化。")
        return

    if G.number_of_nodes() > max_nodes:
        degrees = dict(G.degree())
        nodes_sorted = sorted(degrees, key=degrees.get, reverse=True)[:max_nodes]
        H = G.subgraph(nodes_sorted).copy()
        print(
            f"  节点超过 {max_nodes} 个，仅绘制度数最高的 {max_nodes} 个超像素。"
        )
    else:
        H = G

    communities = [partition[node] for node in H.nodes()]
    n_communities = len(set(communities))
    colors = plt.cm.tab20(np.linspace(0, 1, max(n_communities, 1)))

    plt.figure(figsize=(14, 10))
    pos = nx.spring_layout(H, seed=42, k=0.3)
    node_colors = [colors[partition[node] % len(colors)] for node in H.nodes()]

    nx.draw_networkx_nodes(
        H,
        pos,
        node_color=node_colors,
        node_size=120,
        linewidths=0.5,
        edgecolors="black",
        alpha=0.85,
    )
    nx.draw_networkx_edges(H, pos, alpha=0.2, width=0.5)
    plt.title("Superpixel Similarity Network (Leiden)", fontsize=16, fontweight="bold")
    plt.axis("off")
    plt.tight_layout()

    network_path = output_dir / "leiden_network.png"
    plt.savefig(network_path, dpi=300)
    plt.close()

    print(f"  网络图已保存: {network_path}")


def main():
    feature_path = Path(
        r" "
    )
    output_dir = Path(
        r" "
    )

    print("=" * 80)
    print("Leiden Bootstrap 共识聚类 - 超像素级别")
    print("=" * 80)
    print(f"特征文件: {feature_path}")
    print(f"输出目录: {output_dir}")

    df_features = load_superpixel_features(feature_path)
    df_selected, selection_summary = select_balanced_superpixels(df_features)
    metadata, features_scaled, feature_cols = prepare_superpixel_matrix(df_selected)

    consensus_result = bootstrap_consensus_leiden(
        features_scaled,
        metadata["superpixel_id"].tolist(),
        n_bootstrap=300,
        sample_ratio=0.8,
        k_neighbors=20,
        resolution=1.0,
        k_range=(20, 25),
        resolution_range=(0.6, 1.2),
        random_state=42,
    )

    superpixel_ids = metadata["superpixel_id"].tolist()
    consensus_matrix = consensus_result["consensus"]

    adaptive_min_weight = determine_consensus_threshold(
        consensus_matrix,
        base_threshold=0.2,
        quantile=0.7,
    )

    print("\n" + "=" * 80)
    print("构建共识图参数")
    print("=" * 80)
    print(f"  自适应最小边权阈值: {adaptive_min_weight:.3f}")
    print("  k-nearest neighbours: 40")

    consensus_graph = build_similarity_graph(
        consensus_matrix,
        superpixel_ids,
        k_neighbors=40,
        min_weight=adaptive_min_weight,
        k_cross=None,
        use_jaccard=False,
    )

    target_range = (3, 10)
    resolution_grid = np.linspace(0.2, 1.6, 15)
    (
        final_partition,
        final_modularity,
        final_resolution,
        final_n_communities,
        resolution_table,
    ) = search_leiden_resolutions(
        consensus_graph,
        resolutions=resolution_grid,
        target_range=target_range,
        random_state=42,
    )

    resolution_search_path = output_dir / "leiden_resolution_search.csv"
    if not resolution_table.empty:
        resolution_table.to_csv(resolution_search_path, index=False)
        print("\n" + "=" * 80)
        print("Leiden 分辨率搜索结果已保存")
        print("=" * 80)
        print(f"  {resolution_search_path}")

    if final_partition is None or final_modularity is None:
        raise RuntimeError("共识图为空或未能找到合适的Leiden分辨率。")

    if final_n_communities is None:
        final_n_communities = len(set(final_partition.values()))

    merged = False
    if final_n_communities > target_range[1]:
        print(
            f"社区数量 {final_n_communities} 超出目标上限 {target_range[1]}，"
            "执行层次聚合以合并社区..."
        )
        final_partition, merged = merge_partition_to_target(
            final_partition,
            superpixel_ids,
            consensus_matrix,
            target_k=target_range[1],
        )
        if merged:
            final_n_communities = len(set(final_partition.values()))
            final_modularity = compute_partition_modularity(
                consensus_graph,
                final_partition,
            )
            print(
                f"  合并后社区数: {final_n_communities}, "
                f"模块度: {final_modularity:.4f}"
            )
        else:
            print("  层次聚合未执行（社区数量已不超过目标或数据不足）。")

    if not (target_range[0] <= final_n_communities <= target_range[1]):
        print(
            f"WARNING: 仍未在目标范围 {target_range} 内获得社区数量，"
            f"当前社区数: {final_n_communities}"
        )

    save_consensus_results(
        metadata,
        final_partition,
        consensus_matrix,
        consensus_result["sample_counts"],
        output_dir,
    )
    visualize_consensus_heatmap(
        consensus_matrix,
        metadata,
        final_partition,
        output_dir,
    )
    visualize_network(
        consensus_graph,
        final_partition,
        output_dir,
    )

    print("\n" + "=" * 80)
    print("分析总结")
    print("=" * 80)
    print(f"选中超像素: {len(metadata)}")
    print(f"最终社区数量: {final_n_communities}")
    if final_resolution is not None:
        print(f"最终 Leiden 分辨率 γ: {final_resolution:.3f}")
    print(f"最终模块度: {final_modularity:.4f}")


if __name__ == "__main__":
    main()
