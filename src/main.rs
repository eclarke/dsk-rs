extern crate kmers;
extern crate dsk;
extern crate env_logger;
#[macro_use] extern crate clap;
#[macro_use] extern crate log;
#[macro_use] extern crate error_chain;

use clap::{App, Arg, ArgMatches};

error_chain! {
    links {
        App(::dsk::errors::Error, ::dsk::errors::ErrorKind);
    }
    foreign_links {
        Logging(::log::SetLoggerError);
    }
}

fn main() {
    env_logger::init().expect("Error initializing logger");
    if let Err(ref e) = run() {
        use std::io::Write;
        let stderr = &mut ::std::io::stderr();
        let errmsg = "Error writing to stderr";

        writeln!(stderr, "error: {}", e).expect(errmsg);

        for e in e.iter().skip(1) {
            writeln!(stderr, "caused by: {}", e).expect(errmsg);
        }

        if let Some(backtrace) = e.backtrace() {
            writeln!(stderr, "backtrace: {:?}", backtrace).expect(errmsg);
        }

        ::std::process::exit(1);
    }
}

fn run() -> Result<()> {
    let yaml = load_yaml!("cli.yml");
    let args = App::from_yaml(yaml).get_matches();
    let app = dsk::App::new(args)?;
    info!("\
Starting DSK with parameters: 
\tk:          {}
\tinput:      {}
\tinput fmt:  {:?}
\titerations: {}
\tpartitions: {}
\tworkspace:  {:?}",
    app.k, app.input, app.format, app.iters, app.parts, app.workspace.path()
);

    let max_k = kmers::max_small_k(app.alphabet());
    if app.k <= max_k {
        info!("Using small kmer counter");
        let counter = app.write_small_kmers()?;
        app.count_kmers(&counter)?;
    } else {
        info!("Using large kmer counter");
        let counter = app.write_large_kmers()?;
        app.count_kmers(&counter)?;
    }

    Ok(())
}

