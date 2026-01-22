use ansi_term::Colour;
use csv::ReaderBuilder;
use std::collections::HashMap;
//引入Arc
use std::process;
use std::sync::Arc;

mod arg;
use arg::get_arg;
use std::thread;

fn main() {
    let start_time = std::time::Instant::now();

    let matches = get_arg();

    // 读取参数
    let input = matches.get_one::<String>("input").unwrap();
    let output = matches.get_one::<String>("output").unwrap();
    let threads = matches.get_one::<String>("threads").unwrap();
    let input_type = matches.get_one::<String>("type").unwrap();

    // 校验input_type（仅kallisto/salmon）
    if input_type != "kallisto" && input_type != "salmon" {
        eprintln!(
            "{}",
            Colour::Red.paint("Error: input_type can only be 'kallisto' or 'salmon'!")
        );
        process::exit(1); // 优雅退出，替代panic
    }

    // 新增接收est_counts_matrix_hash
    let (
        sample_names,
        count_matrix_hash,
        tpm_matrix_hash,
        fpkm_matrix_hash,
        est_counts_matrix_hash,
    ) = read_file_2_vec(input, input_type);

    if threads == "1" {
        // 写出count matrix
        write_matrix(&output, count_matrix_hash, "count", &sample_names);
        // 写出tpm matrix
        write_matrix(&output, tpm_matrix_hash, "tpm", &sample_names);
        // 写出fpkm matrix
        write_matrix(&output, fpkm_matrix_hash, "fpkm", &sample_names);
        // 新增：写出est_counts matrix
        write_matrix(&output, est_counts_matrix_hash, "est_counts", &sample_names);
    } else {
        // convert threads 20231007
        let output_arc = Arc::new(output.to_string());
        let output_arc_1 = output_arc.clone();
        let output_arc_2 = output_arc.clone();
        let output_arc_3 = output_arc.clone();
        let output_arc_4 = output_arc.clone(); // 新增：est_counts的output克隆
        let sample_names_arc = Arc::new(sample_names);
        let sample_names_arc_1 = sample_names_arc.clone();
        let sample_names_arc_2 = sample_names_arc.clone();
        let sample_names_arc_3 = sample_names_arc.clone();
        let sample_names_arc_4 = sample_names_arc.clone(); // 新增：est_counts的sample name克隆

        let handle_write_count_matrix = thread::spawn(move || {
            write_matrix(
                &output_arc_1,
                count_matrix_hash,
                "count",
                &sample_names_arc_1,
            )
        });
        let handle_write_tpm_matrix = thread::spawn(move || {
            write_matrix(
                &output_arc_2.clone(),
                tpm_matrix_hash,
                "tpm",
                &sample_names_arc_2,
            )
        });
        let handle_write_fpkm_matrix = thread::spawn(move || {
            write_matrix(&output_arc_3, fpkm_matrix_hash, "fpkm", &sample_names_arc_3)
        });
        // 新增：est_counts矩阵的写入线程
        let handle_write_est_counts_matrix = thread::spawn(move || {
            write_matrix(
                &output_arc_4,
                est_counts_matrix_hash,
                "est_counts",
                &sample_names_arc_4,
            )
        });

        handle_write_count_matrix.join().unwrap();
        handle_write_tpm_matrix.join().unwrap();
        handle_write_fpkm_matrix.join().unwrap();
        handle_write_est_counts_matrix.join().unwrap(); // 新增：等待est_counts线程完成
    }

    println!(
        "{}{}",
        Colour::Green.paint("Done!"),
        Colour::Red.paint("Goodbye!")
    );
    println!(
        "{}{}",
        Colour::Green.paint("Total elapsed time: "),
        Colour::Red.paint(format!("{:?}", start_time.elapsed()))
    );
}

fn read_kallisto(input: &str) -> HashMap<String, Vec<f64>> {
    println!("read kallisto output file: {}", Colour::Green.paint(input));
    let mut map = HashMap::new();
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(input)
        .unwrap();
    let mut total_reads = 0.0;
    for result in reader.records() {
        let record = result.unwrap();
        let gene_id = record.get(0).unwrap().to_string();
        let tpm = record.get(4).unwrap().parse::<f64>().unwrap();
        let eff_counts = record.get(3).unwrap().parse::<f64>().unwrap(); // est_counts即eff_counts
                                                                         // eff_counts = count * (length / eff_length)
                                                                         // 根据公式计算count
                                                                         // count = eff_counts / (length / eff_length)
        let count = eff_counts
            / (record.get(1).unwrap().parse::<f64>().unwrap()
                / record.get(2).unwrap().parse::<f64>().unwrap());
        total_reads += count;
        let mut vec = Vec::new();
        vec.push(tpm); // 索引0: tpm
        vec.push(count); // 索引1: count
        vec.push(eff_counts); // 索引2: est_counts (新增)
        map.insert(gene_id, vec);
    }
    // 根据公式计算fpkm
    // FPKM: Fragments Per Kilobase of exon model per Million mapped fragments
    // FPKM = 10^9 * C / (N * L)
    // C: count number of reads mapped to a gene
    // N: total mapped reads in the experiment
    // L: exonic length in base pairs for a gene
    // let fpkm = 1000000000.0 * count / (total_reads * length);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(input)
        .unwrap();
    for result in reader.records() {
        let record = result.unwrap();
        let gene_id = record.get(0).unwrap().to_string();
        let eff_counts = record.get(3).unwrap().parse::<f64>().unwrap();
        if eff_counts == 0.0 {
            let mut vec = map.get(&gene_id).unwrap().clone();
            vec.push(0.0); // 索引3: fpkm (0值)
            map.insert(gene_id, vec);
            continue;
        }
        let count = eff_counts
            / (record.get(1).unwrap().parse::<f64>().unwrap()
                / record.get(2).unwrap().parse::<f64>().unwrap());
        let length = record.get(1).unwrap().parse::<f64>().unwrap();
        let fpkm = 1000000000.0 * count / (total_reads * length);
        let mut vec = map.get(&gene_id).unwrap().clone();
        vec.push(fpkm); // 索引3: fpkm
        map.insert(gene_id, vec);
    }
    map
}

fn read_salmon(input: &str) -> HashMap<String, Vec<f64>> {
    // Name    Length  EffectiveLength TPM     NumReads
    // A.satnudSFS1A01G007012.1        1413    1279.982        1.913018        41.482
    // A.satnudSFS1A01G007011.1        1503    1280.518        3.720414        80.707
    println!("read salmon output file: {}", Colour::Green.paint(input));
    let mut map = HashMap::new();
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(input)
        .unwrap();
    let mut total_reads = 0.0;
    for result in reader.records() {
        let record = result.unwrap();
        let gene_id = record.get(0).unwrap().to_string();
        let tpm = record.get(3).unwrap().parse::<f64>().unwrap();
        let counts = record.get(4).unwrap().parse::<f64>().unwrap();
        total_reads += counts;
        let mut vec = Vec::new();
        vec.push(tpm); // 索引0: tpm
        vec.push(counts); // 索引1: count
        vec.push(0.0); // 索引2: est_counts (salmon无，填充0，保证格式统一)
        map.insert(gene_id, vec);
    }
    // 根据公式计算fpkm
    // FPKM: Fragments Per Kilobase of exon model per Million mapped fragments
    // FPKM = 10^9 * C / (N * L)
    // C: count number of reads mapped to a gene
    // N: total mapped reads in the experiment
    // L: exonic length in base pairs for a gene
    // let fpkm = 1000000000.0 * count / (total_reads * length);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(input)
        .unwrap();
    for result in reader.records() {
        let record = result.unwrap();
        let gene_id = record.get(0).unwrap().to_string();
        let counts = record.get(4).unwrap().parse::<f64>().unwrap();
        if counts == 0.0 {
            let mut vec = map.get(&gene_id).unwrap().clone();
            vec.push(0.0); // 索引3: fpkm (0值)
            map.insert(gene_id, vec);
            continue;
        }
        let length = record.get(1).unwrap().parse::<f64>().unwrap();
        let fpkm = 1000000000.0 * counts / (total_reads * length);
        let mut vec = map.get(&gene_id).unwrap().clone();
        vec.push(fpkm); // 索引3: fpkm
        map.insert(gene_id, vec);
    }
    map
}

fn write_matrix(
    name: &str,
    out_hash: HashMap<String, Vec<f64>>,
    prefix: &str,
    sample_names: &Vec<String>,
) {
    println!(
        "write {} matrix file: {}",
        Colour::Green.paint(prefix),
        Colour::Green.paint(format!("{}_{}_{}.txt", name, prefix, "matrix"))
    );
    // 优化使用csv库写入
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(format!("{}_{}_{}.txt", name, prefix, "matrix"))
        .unwrap();
    let mut header = Vec::new();
    header.push("gene_id");
    for sample_name in sample_names.iter() {
        header.push(sample_name);
    }
    wtr.write_record(&header).unwrap();
    // 对hashmap进行排序
    let mut out_hash: Vec<_> = out_hash.into_iter().collect();
    out_hash.sort_by(|a, b| a.0.cmp(&b.0));
    // 写出文件
    for (gene_id, counts) in out_hash.iter() {
        let mut line = Vec::new();
        line.push(gene_id.clone());
        counts.iter().for_each(|count| line.push(count.to_string()));
        // 检查行列数匹配，避免错误
        if &line.len() != &header.len() {
            println!(
                "警告: 基因{}的数值列数({})与样本数({})不匹配，跳过该基因",
                gene_id,
                counts.len(),
                sample_names.len()
            );
            continue;
        } else {
            wtr.write_record(&line).unwrap();
        }
    }
    wtr.flush().unwrap();

    println!(
        "write {} matrix file: {} done!",
        Colour::Green.paint(prefix),
        Colour::Green.paint(format!("{}_{}_{}.txt", name, prefix, "matrix"))
    );
}

// 返回多个结果 使用元组（新增est_counts_matrix_hash）
fn read_file_2_vec(
    input: &str,
    input_type: &str,
) -> (
    Vec<String>,
    HashMap<String, Vec<f64>>,
    HashMap<String, Vec<f64>>,
    HashMap<String, Vec<f64>>,
    HashMap<String, Vec<f64>>, // 新增：est_counts矩阵hash
) {
    // input: sample csv
    // ./BB313-01T0001_sfs/abundance.tsv,01T0001_sfs
    // ./BB313-01T0002_sfs/abundance.tsv,01T0002_sfs
    let mut sample_names = Vec::new();
    let mut count_matrix_hash = HashMap::new();
    let mut tpm_matrix_hash = HashMap::new();
    let mut fpkm_matrix_hash = HashMap::new();
    let mut est_counts_matrix_hash = HashMap::new(); // 新增：est_counts矩阵hash

    let mut samples = ReaderBuilder::new()
        .delimiter(',' as u8)
        .has_headers(false)
        .from_path(input)
        .unwrap();
    // 读取每个样本的kallisto/salmon输出文件
    for sample in samples.records() {
        let record = sample.unwrap();
        let sample_name = record.get(1).unwrap().to_string();
        sample_names.push(sample_name);
        let file = record.get(0).unwrap().to_string();
        let map = if input_type == "salmon" {
            read_salmon(&file)
        } else {
            read_kallisto(&file)
        };
        // 优化代码 2023-05-08
        for gene_id in map.keys() {
            // count (索引1)
            let vec = count_matrix_hash
                .entry(gene_id.to_string())
                .or_insert(vec![]);
            vec.push(map.get(gene_id).expect("读取count失败")[1]);

            // tpm (索引0)
            let vec = tpm_matrix_hash.entry(gene_id.to_string()).or_insert(vec![]);
            vec.push(map.get(gene_id).expect("读取tpm失败")[0]);

            // fpkm (索引3)
            let vec = fpkm_matrix_hash
                .entry(gene_id.to_string())
                .or_insert(vec![]);
            vec.push(map.get(gene_id).expect("读取fpkm失败")[3]);

            // est_counts (索引2，新增)
            let vec = est_counts_matrix_hash
                .entry(gene_id.to_string())
                .or_insert(vec![]);
            vec.push(map.get(gene_id).expect("读取est_counts失败")[2]);
        }
    }
    (
        sample_names,
        count_matrix_hash,
        tpm_matrix_hash,
        fpkm_matrix_hash,
        est_counts_matrix_hash, // 新增返回值
    )
}
