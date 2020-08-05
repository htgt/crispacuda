#include "crispacuda.h"
#include "mongoose.h"
#include "options.h"
#include "ots.h"
#include "seq.h"
#include <sstream>
#include <vector>
using namespace std;

void parse_request(struct mg_str req, vector<uint64_t> &ids) {
    char *query, *id;
    query = (char*)malloc(req.len + 1);
    memset(query, 0, req.len + 1);
    strncpy(query, req.p, req.len);
    while ( (id = strsep(&query, "\n")) != NULL ) {
        ids.push_back(atol(id));
    }
}

static void handle_request(struct mg_connection *c, int ev, void *p) {
    userdata_t *data = (struct userdata_t*)c->user_data;
    if ( ev == MG_EV_HTTP_REQUEST ) {
        struct http_message *hm = (struct http_message *)p;
        
        if( mg_vcmp(&hm->uri, "/search") == 0 ) {
            vector<uint64_t> ids;
            parse_request(hm->body, ids);
            stringstream results;
            for(uint64_t id : ids) {
                crispr_t crispr;
                if ( id <= data->metadata.offset || id > data->metadata.offset + data->metadata.num_seqs ) {
                    continue;
                }
                crispr.id = id;
                crispr.seq = data->h_crisprs[id - data->metadata.offset - 1];
                crispr.rev_seq = revcom(crispr.seq, data->metadata.seq_length);
                calc_off_targets(results, data, crispr);
            }
            const string tmp = results.str();
            struct mg_str response = mg_mk_str(tmp.c_str());
            mg_send_head(c, 200, response.len, "Content-Type: text/plain");
            mg_printf(c, "%.*s", (int)response.len, response.p);
        } else if( mg_vcmp(&hm->uri, "/favicon.ico") == 0 ) {
            mg_http_serve_file(c, hm, "favicon.ico", mg_mk_str("image/ico"), mg_mk_str(""));
        } else {
            mg_http_serve_file(c, hm, "index.htm", mg_mk_str("text/html"), mg_mk_str(""));
        }
    }
}

void run_server(char *port, userdata_t userdata) {
    struct mg_mgr mgr;
    struct mg_connection *c;
    mg_mgr_init(&mgr, NULL);
    
    c = mg_bind(&mgr, port, handle_request);
    c->user_data = &userdata;
    mg_set_protocol_http_websocket(c);
    for(;;) {
        mg_mgr_poll(&mgr, 1000000);
    }
    mg_mgr_free(&mgr);
}
