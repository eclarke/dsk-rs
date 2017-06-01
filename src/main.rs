extern crate kmers;
extern crate dsk;
extern crate slog_term;
extern crate slog_async;
#[macro_use] extern crate slog;
#[macro_use] extern crate clap;
#[macro_use] extern crate error_chain;

use slog::Drain;
use clap::App;

error_chain! {
    links {
        App(::dsk::errors::Error, ::dsk::errors::ErrorKind);
    }
}

fn main() {
    // env_logger::init().expect("Error initializing logger");
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

    let decorator = slog_term::TermDecorator::new().build();
    let drain = slog_term::FullFormat::new(decorator).build().fuse();
    let drain = slog_async::Async::new(drain).build().fuse();

    let log = slog::Logger::root(drain, o!());

    let yaml = load_yaml!("cli.yml");
    let args = App::from_yaml(yaml).get_matches();
    let app = dsk::App::new(args, log.new(o!()))?;
    info!(log, "DSK started"; "k"=>&app.k, "input"=>&app.input);
//     info!("\
// Starting DSK with parameters: 
// \tk:          {}
// \tinput:      {}
// \tinput fmt:  {:?}
// \titerations: {}
// \tpartitions: {}
// \tworkspace:  {:?}",
//     app.k, app.input, app.format, app.iters, app.parts, app.workspace.path()
// );

    let max_k = kmers::max_small_k(app.alphabet());
    if app.k <= max_k {
        info!(log, "Writing small kmers to disk");
        let counter = app.write_small_kmers()?;
        info!(log, "Counting kmers");
        let map = app.count_kmers(&counter)?;
        info!(log, "Writing map to disk");
        app.write_map(map)?
    } else {
        info!(log, "Writing large kmers to disk");
        let counter = app.write_large_kmers()?;
        info!(log, "Counting kmers");
        let map = app.count_kmers(&counter)?;
        info!(log, "Writing map to disk");
        app.write_map(map)?
    }
    info!(log, "Finished");

    Ok(())
}

