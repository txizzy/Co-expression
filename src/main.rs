use std::fs::{File, remove_file};
use std::io::Write;
use std::io::{BufRead, BufReader};
use statrs::distribution::{StudentsT,ContinuousCDF};
use polars::prelude::*;
use clap::{App, Arg, ArgMatches};
use std::collections::HashMap;
use rayon::prelude::*;
use std::sync::Mutex;

fn args_input() -> ArgMatches<'static>{
    let args = App::new("Calculate the pearson's correlation by rust")
        .version("0.1.0")
        .author("Xizzy (txizzy#gmail.com)")
        .about("Usage:\nPearson -i file1 -I file2 -o out.xls\n")
        .arg(
            Arg::with_name("infile1")
                .long("infile1")
                .short("i")
                .help("Input the gene expression file1. Use a table to split the columns. The first column should contain the gene names and a header is required.")
                .takes_value(true)
                .required(true))
        .arg(
            Arg::with_name("infile2")
                .long("infile2")
                .short("I")
                .help("Input the gene expression file2, same requirements as file 1.")
                .takes_value(true)
                .required(true)
        )
        .arg(
            Arg::with_name("output")
            .long("output")
            .short("o")
            .help("Output filename")
            .takes_value(true)
            .required(true)
        )
        .arg(
            Arg::with_name("corr")
            .long("corr")
            .short("c")
            .help("Cutoff of correlation's absolute")
            .takes_value(true)
            .default_value("0.95")
        )
        .arg(
            Arg::with_name("pvalue")
            .long("pvalue")
            .short("p")
            .help("Cutoff of pvalue, default is 0.05")
            .takes_value(true)
            .default_value("0.05")
        )
        .arg(
            Arg::with_name("threads")
            .long("threads")
            .short("t")
            .help("Number of threads for parallel computation")
            .takes_value(true)
            .default_value("4")
        )
        .get_matches();
    args
}


fn main() {
    let args = args_input();
    let infile1 = args.value_of("infile1").expect("Can not found the file1");
    let infile2 = args.value_of("infile2").expect("Can not found the file2");
    let output = args.value_of("output").expect("Can not found the output");
    let pvalue: f64 = args.value_of("pvalue").unwrap().parse().unwrap();
    let corr: f64 = args.value_of("corr").unwrap().parse().unwrap();
    let threads: usize = args.value_of("threads").unwrap().parse().unwrap();
    go_polars(infile1, infile2, output, pvalue, corr, threads);
}

fn go_polars(infile1: &str, infile2 :&str, output: &str, pvalue : f64, corr : f64, threads: usize){
    let file1 = File::open(infile1).expect("could not open file1");
    let file2 = File::open(infile2).expect("could not open file2");
    //let mut outfile = File::create(output).expect("could not create output file");

    let read_tsv = |file: File| {
        CsvReader::new(file)
            .infer_schema(None)
            .has_header(true)
            .with_delimiter(b'\t')
            .finish()
            .unwrap()
    };

    let df1 = read_tsv(file1);
    let df2 = read_tsv(file2);

    let nsize = df1.width() - 1 ;
    
    // make gene list for echo dataframe
    let get_key = |ser :&Series | -> Vec<String> {ser.iter().map(|item| 
    item.to_string().replace("\"", "")).collect()};
    let df1_key: Vec<String> = get_key(&df1[0]);

    let df2_key: Vec<String> = get_key(&df2[0]);

    // get expression matrix
    let get_exp = |df : &DataFrame, nsize: usize| -> Vec<Vec<f64>>{
        df.select_by_range(1..nsize+1).unwrap().iter().map(|df_df| df_df.f64().unwrap().into_iter().flatten().collect()).collect()
    };
    let df1_exp: Vec<Vec<f64>> = get_exp(&df1, nsize); 

    //transposed_data 
    let mut df1_exp2 : Vec<Vec<f64>> =  vec![vec![0.0;  df1_exp.len()]; df1_exp[0].len()];
    for (row_idx, row) in df1_exp.iter().enumerate() {
        for (col_idx, &value) in row.iter().enumerate() {
            df1_exp2[col_idx][row_idx] = value;
        }
    }
    let df2_exp: Vec<Vec<f64>> = get_exp(&df2, nsize); 
    
    //transposed_data 
    let mut df2_exp2 : Vec<Vec<f64>> =  vec![vec![0.0;  df2_exp.len()]; df2_exp[0].len()];
    for (row_idx, row) in df2_exp.iter().enumerate() {
        for (col_idx, &value) in row.iter().enumerate() {
            df2_exp2[col_idx][row_idx] = value;
        }
    }

    //chunk df1
    let real_threads = if df1_exp2.len() < (threads * 2) {
        println!("Number of nrows in infile1 are too small, only use 1 thread!");
        1
    } else {
        threads
    };
    let df1_exp_chunks = df1_exp2.chunks(df1_exp2.len() / real_threads + 1);
    let df1_key_chunks = df1_key.chunks(df1_key.len() / real_threads + 1);
    let mut df1_exp_dict: HashMap<usize, Vec<Vec<f64>>> = HashMap::new();
    let mut df1_key_dict: HashMap<usize, Vec<String>> = HashMap::new();
    for (idx,chunk) in df1_exp_chunks.enumerate(){
        df1_exp_dict.insert(idx, chunk.to_vec());
    }
    for (idx, chunk) in df1_key_chunks.enumerate(){
        df1_key_dict.insert(idx, chunk.to_vec());
    }
   
     // 创建指定数量的线程
     (0..real_threads).into_par_iter().for_each(|i| {
        let outname = format!("{}_{}.temp", output,i);
        //println!("{:?}\n{:?}\n{:?}\n{}",df1_exp_dict[&i],df1_key_dict[&i],df2_exp2, outname);
        call_pearson(df1_exp_dict[&i].clone(), df2_exp2.clone(), df1_key_dict[&i].clone(), df2_key.clone(), corr, pvalue, outname);
    });


    let output_mutex = Mutex::new(File::create(output).unwrap());
    (0..real_threads).into_par_iter().for_each(|id| {
        let outname = format!("{}_{}.temp", output, id);
        let file = File::open(&outname).unwrap();
        let reader = BufReader::new(file);
        let buffer_size = 100000;
        let mut lines = Vec::with_capacity(buffer_size);

        for line in reader.lines() {
            lines.push(line.unwrap() + "\n");

            if lines.len() >= buffer_size {
                let mut outfile = output_mutex.lock().unwrap();
                outfile.write_all(lines.join("").as_bytes()).unwrap();
                lines.clear();
            }
        }

        if !lines.is_empty() {
            let mut outfile = output_mutex.lock().unwrap();
            outfile.write_all(lines.join("").as_bytes()).unwrap();
        }

        remove_file(outname).unwrap();
    });
    
}


fn call_pearson(df1_exp2: Vec<Vec<f64>>, df2_exp2: Vec<Vec<f64>>, 
                df1_key: Vec<String>, df2_key: Vec<String> ,corr: f64, pvalue: f64, output: String ) {
    let mut outfile = File::create(output).expect("could not create output file");
    for x in 0..df1_exp2.len(){
        for y in 0..df2_exp2.len(){
            let (correlations, p) = pearson_correlation( &df1_exp2[x], &df2_exp2[y] );
            if (p < pvalue) && (correlations.abs() > corr){
                let result = format!("{}\t{}\t{}\t{}\n", df1_key[x], df2_key[y], correlations, p);
                let _ = outfile.write(result.as_bytes());
            }
        }
    }
}

// Calculate Pearson
fn pearson_correlation(x: &[f64], y: &[f64]) -> (f64,f64) {
    let n = x.len();
    let sum_x: f64 = x.iter().sum();
    let sum_y: f64 = y.iter().sum();
    if sum_x > 0.0 && sum_y > 0.0 {
        let sum_xy: f64 = x.iter().zip(y.iter()).map(|(&xi, &yi)| xi * yi).sum();
        let sum_x_squared: f64 = x.iter().map(|&xi| xi * xi).sum();
        let sum_y_squared: f64 = y.iter().map(|&yi| yi * yi).sum();

        let numerator = n as f64 * sum_xy - sum_x * sum_y;
        let denominator = ((n as f64 * sum_x_squared - sum_x * sum_x) * (n as f64 * sum_y_squared - sum_y * sum_y)).sqrt();
        let correlations = if numerator / denominator > 1.0 {
            1.0 as f64
        }else{
            numerator / denominator 
        };
        let t = correlations.abs() / ((1.0 as f64 - correlations*correlations)/ (n - 2)  as f64).sqrt();
        let student_t = StudentsT::new(0.0, 1.0, (n - 2) as f64).unwrap();
        let p = (1.0 -student_t.cdf(t.abs())) * 2.0;
        (correlations, p)
    }else{
        (0.0, 1.0)
    }
}
