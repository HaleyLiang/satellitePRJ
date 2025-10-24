#include <iostream>
#include "src/Funs.h"
#include "cpp-httplib/httplib.h"
#include "src/visibility_calculator.h"
#include "src/visibility_service.h"
#include "src/delay_service.h"
#include "src/link_margin_service.h"
#include "src/communication_quality_service.h"
#include "src/interference_assessment_service.h"
#include "src/link_assessment_service.h"
#include "src/simplified_setup_time_service.h"

int main() {
    // 实例化服务器
    httplib::Server svr;

    // 设置一个路由
    svr.Get("/", [](const httplib::Request&, httplib::Response& res) {
        res.set_content("Hello World from MinGW!", "text/plain");
    });

    // 另一个显示信息的路由
    svr.Get("/stop", [&](const httplib::Request&, httplib::Response& res) {
        res.set_content("Stopping the server...", "text/plain");
        svr.stop(); // 通知服务器停止
    });

    // 可见性计算接口(多个卫星单个地面站)
    svr.Post("/check_visibility", [](const httplib::Request& req, httplib::Response& res) {
        try {
            auto jsonBody = json::parse(req.body);
            auto results = calculateVisibility(jsonBody["terminal"], jsonBody["satellites"]);
            res.set_content(results.dump(), "application/json");
        } catch (const std::exception& e) {
            res.status = 400;
            res.set_content(json{{"error", e.what()}}.dump(), "application/json");
        }
    });

    // 可见性计算接口（多个卫星和多个地面站）
    svr.Post("/check_multi_visibility", [](const httplib::Request& req, httplib::Response& res) {
        try {
            auto jsonBody = json::parse(req.body);

            // 验证请求格式
            if (!jsonBody.contains("terminals") || !jsonBody.contains("satellites")) {
                res.status = 400;
                res.set_content(json{{"error", "Missing required fields: terminals or satellites"}}.dump(), "application/json");
                return;
            }

            auto results = VisibilityService::calculateMultiVisibility(
                jsonBody["terminals"],
                jsonBody["satellites"]
            );

            res.set_content(results.dump(), "application/json");
        } catch (const std::exception& e) {
            res.status = 400;
            res.set_content(json{{"error", e.what()}}.dump(), "application/json");
        }
    });

    // 时延计算接口
    svr.Post("/calculate_delays", [](const httplib::Request& req, httplib::Response& res) {
        try {
            auto jsonBody = json::parse(req.body);

            // 验证请求格式
            if (!jsonBody.contains("terminals") || !jsonBody.contains("satellites")) {
                res.status = 400;
                res.set_content(json{{"error", "Missing required fields: terminals or satellites"}}.dump(), "application/json");
                return;
            }

            auto results = DelayService::calculateMultiDelay(
                jsonBody["terminals"],
                jsonBody["satellites"]
            );

            res.set_content(results.dump(), "application/json");
        } catch (const std::exception& e) {
            res.status = 400;
            res.set_content(json{{"error", e.what()}}.dump(), "application/json");
        }
    });

    // 链路余量计算接口
    svr.Post("/calculate_link_margins", [](const httplib::Request& req, httplib::Response& res) {
        try {
            auto jsonBody = json::parse(req.body);

            // 验证请求格式
            if (!jsonBody.contains("links")) {
                res.status = 400;
                res.set_content(json{{"error", "Missing required field: links"}}.dump(), "application/json");
                return;
            }

            auto results = LinkMarginService::calculateAllLinkMargins(jsonBody["links"]);
            res.set_content(results.dump(), "application/json");
        } catch (const std::exception& e) {
            res.status = 400;
            res.set_content(json{{"error", e.what()}}.dump(), "application/json");
        }
    });

    // 通信质量计算接口
    svr.Post("/calculate_communication_quality", [](const httplib::Request& req, httplib::Response& res) {
        try {
            auto jsonBody = json::parse(req.body);

            // 验证请求格式
            if (!jsonBody.contains("links")) {
                res.status = 400;
                res.set_content(json{{"error", "Missing required field: links"}}.dump(), "application/json");
                return;
            }

            auto results = CommunicationQualityService::calculateAllLinkQualities(jsonBody["links"]);
            res.set_content(results.dump(), "application/json");
        } catch (const std::exception& e) {
            res.status = 400;
            res.set_content(json{{"error", e.what()}}.dump(), "application/json");
        }
    });

    // 受干扰情况评估接口
    svr.Post("/assess_interference", [](const httplib::Request& req, httplib::Response& res) {
        try {
            auto jsonBody = json::parse(req.body);

            // 验证请求格式
            if (!jsonBody.contains("links")) {
                res.status = 400;
                res.set_content(json{{"error", "Missing required field: links"}}.dump(), "application/json");
                return;
            }

            auto results = InterferenceAssessmentService::assessAllLinksInterference(jsonBody["links"]);
            res.set_content(results.dump(), "application/json");
        } catch (const std::exception& e) {
            res.status = 400;
            res.set_content(json{{"error", e.what()}}.dump(), "application/json");
        }
    });

    // 链路评估接口
    svr.Post("/assess_links", [](const httplib::Request& req, httplib::Response& res) {
        try {
            auto jsonBody = json::parse(req.body);

            // 验证请求格式
            if (!jsonBody.contains("links")) {
                res.status = 400;
                res.set_content(json{{"error", "Missing required field: links"}}.dump(), "application/json");
                return;
            }

            auto results = LinkAssessmentService::assessAllLinks(jsonBody["links"]);
            res.set_content(results.dump(), "application/json");
        } catch (const std::exception& e) {
            res.status = 400;
            res.set_content(json{{"error", e.what()}}.dump(), "application/json");
        }
    });

    // 建链时长计算接口
    svr.Post("/calculate_setup_times", [](const httplib::Request& req, httplib::Response& res) {
        try {
            auto jsonBody = json::parse(req.body);

            // 验证请求格式
            if (!jsonBody.contains("links")) {
                res.status = 400;
                res.set_content(json{{"error", "Missing required field: links"}}.dump(), "application/json");
                return;
            }

            auto results = SimplifiedSetupTimeService::calculateAllSetupTimes(jsonBody["links"]);
            res.set_content(results.dump(), "application/json");
        } catch (const std::exception& e) {
            res.status = 400;
            res.set_content(json{{"error", e.what()}}.dump(), "application/json");
        }
    });

    std::cout << "Server started at http://localhost:8080" << std::endl;
    std::cout << "Press Ctrl+C to exit, or visit /stop" << std::endl;

    // 启动服务器，监听所有地址的 8080 端口
    // listen() 会阻塞，直到服务器被停止
    svr.listen("0.0.0.0", 8080);

    std::cout << "Server has stopped." << std::endl;
    return 0;
}