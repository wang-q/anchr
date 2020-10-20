use clap::*;
use std::collections::HashMap;
use tera::{Context, Tera};

// Create clap subcommand arguments
pub fn make_subcommand<'a, 'b>() -> App<'a, 'b> {
    SubCommand::with_name("template")
        .about("Creates Bash scripts")
        .after_help(
            r#"
* Bash scripts:

    * `select_col.sh`
    * `result_stat.sh`
    * `result_dist.sh`
    * `BS.sh`
    * `bootstrap.sh`
    * `filter_result.sh`

"#,
        )
        .arg(Arg::with_name("all").long("all").help("Create all scripts"))
        .arg(
            Arg::with_name("threshold")
                .long("threshold")
                .takes_value(true)
                .default_value("1825")
                .empty_values(false)
                .help("Threshold of time to be censored"),
        )
        .arg(
            Arg::with_name("parallel")
                .long("parallel")
                .takes_value(true)
                .default_value("8")
                .empty_values(false)
                .help("Number of threads"),
        )
        .arg(
            Arg::with_name("sample")
                .long("sample")
                .takes_value(true)
                .default_value("150")
                .empty_values(false)
                .help("Number of sampled cases for bootstrap"),
        )
        .arg(
            Arg::with_name("count")
                .long("count")
                .takes_value(true)
                .default_value("100")
                .empty_values(false)
                .help("Number of bootstrap tests"),
        )
        .arg(
            Arg::with_name("coxp")
                .long("coxp")
                .takes_value(true)
                .default_value("0.05")
                .empty_values(false)
                .help("P value threshold for Cox PH Model"),
        )
        .arg(
            Arg::with_name("hr")
                .long("hr")
                .takes_value(true)
                .default_value("0.95,1.05")
                .empty_values(false)
                .help("HR threshold for KM estimator"),
        )
        .arg(
            Arg::with_name("kmp")
                .long("kmp")
                .takes_value(true)
                .default_value("1")
                .empty_values(false)
                .help("P value threshold for KM estimator"),
        )
        .arg(
            Arg::with_name("rocauc")
                .long("rocauc")
                .takes_value(true)
                .default_value("0.49,0.51")
                .empty_values(false)
                .help("AUC threshold for ROC analysis"),
        )
        .arg(
            Arg::with_name("rocp")
                .long("rocp")
                .takes_value(true)
                .default_value("1")
                .empty_values(false)
                .help("P value threshold for ROC analysis"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> std::result::Result<(), std::io::Error> {
    // context from args
    let mut opt = HashMap::new();
    opt.insert("threshold", args.value_of("threshold").unwrap());
    opt.insert("parallel", args.value_of("parallel").unwrap());
    opt.insert("sample", args.value_of("sample").unwrap());
    opt.insert("count", args.value_of("count").unwrap());
    opt.insert("coxp", args.value_of("coxp").unwrap());
    opt.insert("hr", args.value_of("hr").unwrap());
    opt.insert("kmp", args.value_of("kmp").unwrap());
    opt.insert("rocauc", args.value_of("rocauc").unwrap());
    opt.insert("rocp", args.value_of("rocp").unwrap());

    let mut context = Context::new();
    context.insert("opt", &opt);

    // create scripts
    if args.is_present("all") {
        gen_univariate(&context)?;

        gen_select_col(&context)?;
    }

    Ok(())
}

fn gen_univariate(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "univariate.R";
    eprintln!("Create {}", outname);

    let template = r###"#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(getopt)
    library(readr)
    library(survival)
    library(pROC)
    library(verification)
})

spec = matrix(
    c(
        "help",
        "h",
        0,
        "logical",
        "brief help message",

        "infile",
        "i",
        1,
        "character",
        "input filename",

        "outfile",
        "o",
        1,
        "character",
        "output filename",

        "valid",
        "v",
        1,
        "double",
        "ratio of valid (non-NA) fields in columns",

        "threshold",
        "t",
        1,
        "integer",
        "threshold of time for censoring",

        "parallel",
        "p",
        1,
        "integer",
        "number of threads"
    ),
    byrow = TRUE,
    ncol = 5
)
opt = getopt(spec)

if (!is.null(opt$help)) {
    cat(getopt(spec, usage = TRUE))
    q(status = 1)
}

if (is.null(opt$infile)) {
    cat("--infile is need\n")
    cat(getopt(spec, usage = TRUE))
    q(status = 1)
}

if (is.null(opt$threshold)) {
    opt$threshold <- {{ opt.threshold }}
}

if (is.null(opt$valid)) {
    opt$valid <- 0.7
}

if (!is.null(opt$outfile)) {
    sink(opt$outfile)
}

# col_type = cols() suppress the output
tbl <- read_tsv(opt$infile, col_type = cols())

# outcome: 0 for bad, 1 for good
# left censored
outcome <-
    ifelse(tbl$time < opt$threshold, ifelse(tbl$status == 1, 0, NA), 1)

# status: 1 for event and 0 for censored
# right censored, only for coxph and K-M plot
tbl$status <- ifelse(tbl$time >= opt$threshold, 0, tbl$status)

do_coxph <- function(tbl, name) {
    if (which(colnames(tbl) == name) <= 3) {
        return()
    }

    if (sum(is.na(tbl[[name]])) / nrow(tbl) > opt$valid) {
        return()
    }

    #----------------------------#
    # coxph
    #----------------------------#
    formula_string = paste('Surv(time, status) ~ ', name, sep = '')
    res_cox <- coxph(as.formula(formula_string), tbl)
    summary <- summary(res_cox)
    coef <- summary$coefficients[1]
    p_value <- summary$coefficients[5]

    #----------------------------#
    # K-M
    #----------------------------#
    predictor <- coef * tbl[[name]]
    score_median <- median(predictor, na.rm = TRUE)
    group <- ifelse(predictor < score_median, 1, 0)

    diff <- tryCatch (
        {
            survdiff(Surv(tbl$time, tbl$status) ~ group, rho = 0)
        },
        error=function(cond) {
            message(cond, "\n")
            return(NA)
        }
    )

    hr <- tryCatch (
        {
            (diff$obs[1]/diff$exp[1]) / (diff$obs[2]/diff$exp[2])
        },
        error=function(cond) {
            message(cond, "\n")
            return(NA)
        }
    )

    kmp <- tryCatch (
        {
            pchisq(diff$chisq, length(diff$n) - 1, lower.tail = F)
        },
        error=function(cond) {
            message(cond, "\n")
            return(NA)
        }
    )

    #----------------------------#
    # ROC
    #----------------------------#
    tbl_tmp <- data.frame(outcome, predictor)
    tbl_roc <-
        tbl_tmp[which(!is.na(tbl_tmp$outcome) &
                         !is.na(tbl_tmp$predictor)), ]

    rocauc <- tryCatch (
        {
            pROC::roc(tbl_roc$outcome,
                tbl_roc$predictor,
                levels = c(0, 1),
                direction = ifelse(hr > 1, ">", "<"),
                na.rm = TRUE)$auc
        },
        error=function(cond) {
            message(cond, "\n")
            return(NA)
        }
    )

    rocp <- tryCatch (
        {
            roc_area <- roc.area(tbl_roc$outcome, group)
            signif(roc_area$p.value, digits=2)
        },
        error=function(cond) {
            message(cond, "\n")
            return(NA)
        }
    )

    sprintf(
        "%s\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n",
        name,
        coef,
        p_value,
        score_median,
        hr,
        kmp,
        rocauc,
        rocp
    )
}

cat("#marker\tcoef\tp\tmedian\thr\tkmp\trocauc\trocp\n")

if (!is.null(opt$parallel)) {
    suppressPackageStartupMessages({
        library("foreach")
        library("doParallel")
        registerDoParallel(cores = opt$parallel)
    })

    foreach(name = colnames(tbl)) %dopar% {
        cat(do_coxph(tbl, name))
    }
    quit() # suppress return values of foreach
} else {
    for (name in colnames(tbl)) {
        cat(do_coxph(tbl, name))
    }
}
"###;

    let rendered = Tera::one_off(template, &context, false).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_select_col(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "select_col.sh";
    eprintln!("Create {}", outname);

    let template = r###"#!/usr/bin/env bash

#----------------------------#
# helper functions
#----------------------------#
set +e

signaled () {
    log_warn Interrupted
    exit 1
}
trap signaled TERM QUIT INT

#----------------------------#
# getopts
#----------------------------#
usage () {
    echo "bash select_col.sh [-f 1-3] <infile> [select]" 1>&2;
    echo "bash select_col.sh -f 1,5 tests/SKCM/simple.data.tsv | datamash check" 1>&2;
    echo "# 6 lines, 2 fields" 1>&2;
    echo "bash select_col.sh -f 1-3 tests/SKCM/meth.tsv.gz tests/SKCM/ucox.05.result.tsv | datamash check" 1>&2;
    echo "# 303 lines, 78 fields" 1>&2;
    exit 1;
}
[ $# -eq 0 ] && usage

while getopts ":f:" opt; do
    case ${opt} in
        f )
            opt_field=$OPTARG
            ;;
        \? )
            echo "Invalid option: $OPTARG" 1>&2
            usage
            ;;
        : )
            echo "Invalid option: $OPTARG requires an argument" 1>&2
            usage
            ;;
    esac
done
shift $((OPTIND - 1))

if [ ! -f "$1" ]; then
    echo "[$1] is not a file" 1>&2; exit;
fi
opt_infile="$1"

if [ "$2" != "" ]; then
    if [ ! -f "$2" ]; then
        echo "[$2] is not a file" 1>&2; exit;
    fi
    opt_select="$2"
fi

#----------------------------#
# run
#----------------------------#
# tmpdir
mytmpdir=$(mktemp -d 2>/dev/null || mktemp -d -t 'mytmpdir')

touch ${mytmpdir}/list

if [ "${opt_field}" != "" ]; then
    echo "${opt_field}" >> ${mytmpdir}/list
fi

if [ "${opt_select}" != "" ]; then
    gzip -dcf "${opt_infile}" |
        head -n 1 |
        tr $'\t' '\n' |
        nl |
        perl -nl -e '
            my @fields = grep {/\w+/} split /\s+/, $_;
            next unless @fields == 2;
            print join qq{\t}, $fields[1], $fields[0];
        ' \
        > ${mytmpdir}/fields.index

    cat "${opt_select}" |
        tsv-select -f 1 \
        > ${mytmpdir}/fields.select

    cat ${mytmpdir}/fields.index |
        grep -Fw -f ${mytmpdir}/fields.select |
        cut -f 2 \
        >> ${mytmpdir}/list
fi

cat ${mytmpdir}/list |
    perl -nl -MAlignDB::IntSpan -e '
        BEGIN {
            our $set = AlignDB::IntSpan->new();
        }

        $set->add($_);

        END {
            my @elements = $set->elements;
            while ( scalar @elements ) {
                my @batching = splice @elements, 0, 5000;

                my $batching_set = AlignDB::IntSpan->new;
                $batching_set->add(@batching);
                print $batching_set->runlist();
            }
        }
    ' \
    > ${mytmpdir}/runlist

count=$(cat ${mytmpdir}/runlist | wc -l)

if [ "$count" -eq "0" ]; then
    echo >&2 "No fields"
elif [ "$count" -eq "1" ]; then
    echo >&2 "Writing fields..."
    gzip -dcf ${opt_infile} |
        tsv-select -f $(cat ${mytmpdir}/runlist)
else
    cat ${mytmpdir}/runlist |
        parallel --no-run-if-empty --line-buffer -k -j 1 --seqreplace ,, "
            echo >&2 'Writing fields ,,...'
            gzip -dcf ${opt_infile} |
                tsv-select -f {} |
                datamash transpose
        " \
        > ${mytmpdir}/selected

    echo >&2 "Writing all fields..."
    cat ${mytmpdir}/selected | datamash transpose
fi

# clean
rm -fr ${mytmpdir}
"###;

    let rendered = Tera::one_off(template, &context, false).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}
